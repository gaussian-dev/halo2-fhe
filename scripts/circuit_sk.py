import os
from random import randint
import argparse
import json
from utils import assign_to_circuit, count_advice_cells_needed_for_poly_range_check, print_advice_cells_info,poly_range_check,CiphertextFormationInput,range_check_ciphertext_formation
from rlwe import rlwe_encryption
from bfv.crt import CRTModuli
from bfv.bfv import BFVCrt
from bfv.discrete_gauss import DiscreteGaussian
from bfv.polynomial import Polynomial


def main(args):
    '''
    ENCRYPTION PHASE - performed outside the circuit.
    '''

    n = args.n
    qis = args.qis
    qis = json.loads(qis)
    t = args.t

    crt_moduli = CRTModuli(qis)
    sigma = 3.2
    discrete_gaussian = DiscreteGaussian(sigma)
    bfv_crt = BFVCrt(crt_moduli, n, t, discrete_gaussian)

    # Perform encryption of m in each CRT basis
    s = bfv_crt.SecretKeyGen()
    e = bfv_crt.bfv_q.rlwe.SampleFromErrorDistribution()
    ais = []
    for i in range(len(crt_moduli.qis)):
        ais.append(bfv_crt.bfv_qis[i].rlwe.Rq.sample_polynomial())

    m = bfv_crt.bfv_q.rlwe.Rt.sample_polynomial()

    ctis = bfv_crt.SecretKeyEncrypt(s, ais, e, m)

    # Sanity check for valid decryption
    message_prime = bfv_crt.Decrypt(s, ctis)

    assert m == message_prime

    # k1 = [QM]t namely the scaled message polynomial
    k1 = m.scalar_mul(crt_moduli.q)
    k1.reduce_coefficients_by_modulus(t)

    # `p` is the modulus of the prime field of the circuit
    p = 21888242871839275222246405745257275088548364400416034343698204186575808495617

    # `r2is` are the polynomials r2i for each i-th CRT basis.
    r2is = []

    # `r1is` are the polynomials r1i for each i-th CRT basis.
    r1is = []

    # `k0is` are the negative multiplicative inverses of t modulo each qi.
    k0is = []

    # `ct0is` are the polynomials ct0i for each CRT basis.
    ct0is = []

    # `ct0is_hat` are the polynomials ct0i_hat for each CRT basis.
    ct0is_hat = []

    '''
    SETUP PHASE - performed outside the circuit
    For each CRT basis, we need to compute the polynomials r1i and r2i (check this doc for more details: https://hackmd.io/@gaussian/r1W98Kqqa)
    '''

    cyclo = [1] + [0] * (n - 1) + [1]
    cyclo = Polynomial(cyclo)

    for i, cti in enumerate(ctis):

        ct0i = cti[0]
        ai = cti[1].scalar_mul(-1)

        #generates the parameters which are used to calculate ct0i
        ct0i_elements = rlwe_encryption(ct0i, qis[i], s, ai, e, k1, cyclo, t, n)

        r2is.append(ct0i_elements[0])
        r1is.append(ct0i_elements[1])
        k0is.append(ct0i_elements[2])
        ct0is.append(ct0i_elements[3])
        ct0is_hat.append(ct0i_elements[4])

    # `r1_bounds` are the bounds for the coefficients of r1i for each CRT basis
    r1_bounds = []

    # `r2_bounds` are the bounds for the coefficients of r2i for each CRT basis
    r2_bounds = []

    # initiate counters for the number of advice cells needed for each constraint phase
    phase_0_assignment_advice_cell_count = 0
    phase_1_range_check_advice_cell_count = 0
    phase_1_eval_at_gamma_constraint_advice_cell_count = 0
    phase_1_encryption_constraint_advice_cell_count = 0

    '''
    CIRCUIT - PHASE 0 - ASSIGNMENT
    '''

    # Every assigned value must be an element of the field Zp. Negative coefficients `-z` are assigned as `p - z`
    s_assigned = assign_to_circuit(s, p)
    e_assigned = assign_to_circuit(e, p)
    k1_assigned = assign_to_circuit(k1, p)

    phase_0_assignment_advice_cell_count += len(s_assigned.coefficients)
    phase_0_assignment_advice_cell_count += len(e_assigned.coefficients)
    phase_0_assignment_advice_cell_count += len(k1_assigned.coefficients)

    r1is_assigned = []
    r2is_assigned = []
    ais_assigned = []
    ct0is_assigned = []

    for i in range(len(ctis)):
        r1i_assigned = assign_to_circuit(r1is[i], p)
        r2i_assigned = assign_to_circuit(r2is[i], p)
        r1is_assigned.append(r1i_assigned)
        r2is_assigned.append(r2i_assigned)

        phase_0_assignment_advice_cell_count += len(r1i_assigned.coefficients)
        phase_0_assignment_advice_cell_count += len(r2i_assigned.coefficients)

        ai_assigned = assign_to_circuit(ais[i], p)
        ct0i_assigned = assign_to_circuit(ct0is[i], p)
        ais_assigned.append(ai_assigned)
        ct0is_assigned.append(ct0i_assigned)

        phase_0_assignment_advice_cell_count += len(ai_assigned.coefficients)
        phase_0_assignment_advice_cell_count += len(ct0i_assigned.coefficients)

    # For the sake of simplicity, we generate a random challenge here
    gamma = randint(0, 1000)

    '''
    CIRCUIT - PHASE 1
    '''

    # Every assigned value must be an element of the field Zp. Negative coefficients `-z` are assigned as `p - z`
    qi_constants = []
    k0i_constants = []

    for i in range(len(ctis)):
        qi_constants.append(qis[i])

        k0i_constant = assign_to_circuit(Polynomial([k0is[i]]), p).coefficients[0]
        k0i_constants.append(k0i_constant)

    '''
    CIRCUIT - PHASE 1 - RANGE CHECK
    '''

    lookup_bits = 8

    # constraint. The coefficients of s should be in the range [-1, 0, 1]
    s_bound = 1
    poly_range_check(s,s_assigned,int(s_bound),p)

    phase_1_range_check_advice_cell_count += count_advice_cells_needed_for_poly_range_check(s_assigned, 2*s_bound + 1, lookup_bits)

    # constraint. The coefficients of e should be in the range [-B, B] where B is the upper bound of the discrete Gaussian distribution
    b = int(discrete_gaussian.z_upper)
    poly_range_check(e,e_assigned,b,p)

    phase_1_range_check_advice_cell_count += count_advice_cells_needed_for_poly_range_check(e_assigned, 2*b + 1, lookup_bits)

    # constraint. The coefficients of k1 should be in the range [-(t-1)/2, (t-1)/2]
    k1_bound = int((t - 1) / 2)
    poly_range_check(k1,k1_assigned,k1_bound,p)

    phase_1_range_check_advice_cell_count += count_advice_cells_needed_for_poly_range_check(k1_assigned, 2*k1_bound + 1, lookup_bits)

    s_at_gamma_assigned = s_assigned.evaluate(gamma)
    phase_1_eval_at_gamma_constraint_advice_cell_count += len(s_assigned.coefficients) * 2 - 1

    e_at_gamma_assigned = e_assigned.evaluate(gamma)
    phase_1_eval_at_gamma_constraint_advice_cell_count += len(e_assigned.coefficients) * 2 - 1

    k1_at_gamma_assigned = k1_assigned.evaluate(gamma)
    phase_1_eval_at_gamma_constraint_advice_cell_count += len(k1_assigned.coefficients) * 2 - 1

    cyclo_at_gamma = cyclo.evaluate(gamma)
    cyclo_at_gamma_assigned = assign_to_circuit(Polynomial([cyclo_at_gamma]), p).coefficients[0]
    n_bits_N = n.bit_length()

    # This corresponds to `load_rlc_cache`
    phase_1_eval_at_gamma_constraint_advice_cell_count += (n_bits_N - 1) * 4
    # This corresponds to `add`
    phase_1_eval_at_gamma_constraint_advice_cell_count += 4

    for i in range(len(ctis)):
        args = CiphertextFormationInput(
            qis=qis[i],
            n=n,
            b=b,
            t=t,
            ctis=ct0is[i],
            ais_pk0=ais[i],
            s_u=s,
            e=e,
            r2is=r2is[i],
            k0is=k0is[i],
            k1=k1,
            ctis_hat=ct0is_hat[i],
            r1is=r1is[i]
        )

        #conducts a range check at each step of forming ct0i
        range_check_ciphertext_formation(args)

        # constraint. The coefficients of r2i should be in the range [-(qi-1)/2, (qi-1)/2]
        r2i_bound = int((qis[i] - 1) / 2)
        r2_bounds.append(r2i_bound)
        poly_range_check(r2is[i],r2is_assigned[i],p,r2i_bound)

        phase_1_range_check_advice_cell_count += count_advice_cells_needed_for_poly_range_check(r2is_assigned[i], 2*r2i_bound + 1, lookup_bits)

        # constraint. The coefficients of (ct0i - ct0i_hat - r2i * cyclo) / qi = r1i should be in the range $[\frac{- ((N+2) \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_{0,i}|)}{q_i}, \frac{(N+2) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_{0,i}|}{q_i}]$
        r1i_bound = (int((qis[i] - 1) / 2) * (n + 2) + b + int((t - 1) / 2) * abs(k0is[i])) / qis[i]
        # round bound to the nearest integer
        r1i_bound = int(r1i_bound)
        r1_bounds.append(r1i_bound)
        poly_range_check(r1is[i],r1is_assigned[i],p,r1i_bound)


        phase_1_range_check_advice_cell_count += count_advice_cells_needed_for_poly_range_check(r1is_assigned[i], 2*r1i_bound + 1, lookup_bits)

        '''
        CIRCUIT - PHASE 1 - EVALUATION AT GAMMA CONSTRAINT
        '''

        r1i_gamma_assigned = r1is_assigned[i].evaluate(gamma)
        r2i_gamma_assigned = r2is_assigned[i].evaluate(gamma)

        phase_1_eval_at_gamma_constraint_advice_cell_count += len(r1is_assigned[i].coefficients) * 2 - 1
        phase_1_eval_at_gamma_constraint_advice_cell_count += len(r2is_assigned[i].coefficients) * 2 - 1

        ai_gamma_assigned = ais_assigned[i].evaluate(gamma)
        ct0i_gamma_assigned = ct0is_assigned[i].evaluate(gamma)

        phase_1_eval_at_gamma_constraint_advice_cell_count += len(ais_assigned[i].coefficients) * 2 - 1
        phase_1_eval_at_gamma_constraint_advice_cell_count += len(ct0is_assigned[i].coefficients) * 2 - 1

        '''
        CIRCUIT - PHASE 1 - CORRECT ENCRYPTION CONSTRAINT
        '''

        lhs = ct0i_gamma_assigned
        rhs = (ai_gamma_assigned * s_at_gamma_assigned + e_at_gamma_assigned + (k1_at_gamma_assigned * k0i_constants[i]) + (r1i_gamma_assigned * qi_constants[i]) + (r2i_gamma_assigned * cyclo_at_gamma_assigned))
        phase_1_encryption_constraint_advice_cell_count += 16

        assert lhs % p == rhs % p

        '''
        VERIFICATION PHASE
        '''

        cyclo_at_gamma = cyclo.evaluate(gamma)
        cyclo_at_gamma_assigned_expected = assign_to_circuit(Polynomial([cyclo_at_gamma]), p).coefficients[0]
        assert cyclo_at_gamma_assigned == cyclo_at_gamma_assigned_expected

        ai_gamma_assigned_expected = ais_assigned[i].evaluate(gamma)
        assert ai_gamma_assigned == ai_gamma_assigned_expected

        ct0i_gamma_assigned_expected = ct0is_assigned[i].evaluate(gamma)
        assert ct0i_gamma_assigned == ct0i_gamma_assigned_expected

        assert qis[i] == qi_constants[i]

        k0i_assigned_expected = assign_to_circuit(Polynomial([k0is[i]]), p).coefficients[0]
        assert k0i_constants[i] == k0i_assigned_expected

    total_advice_cell_count = phase_0_assignment_advice_cell_count + phase_1_range_check_advice_cell_count + phase_1_eval_at_gamma_constraint_advice_cell_count + phase_1_encryption_constraint_advice_cell_count

    print_advice_cells_info(total_advice_cell_count, phase_0_assignment_advice_cell_count, phase_1_range_check_advice_cell_count, phase_1_eval_at_gamma_constraint_advice_cell_count, phase_1_encryption_constraint_advice_cell_count)
    # ais and ct0is need to be parsed such that their coefficients are in the range [0, p - 1]
    # we don't call them assigned because they are never assigned to the circuit
    ais_in_p = [assign_to_circuit(ai, p) for ai in ais]
    ct0is_in_p = [assign_to_circuit(ct0i, p) for ct0i in ct0is]
    
    # Parse the inputs into a JSON format such this can be used as input for the (real) circuit
    json_input = {
        "s": [str(coef) for coef in s_assigned.coefficients],
        "e": [str(coef) for coef in e_assigned.coefficients],
        "k1": [str(coef) for coef in k1_assigned.coefficients],
        "r2is": [[str(coef) for coef in r2i.coefficients] for r2i in r2is_assigned],
        "r1is": [[str(coef) for coef in r1i.coefficients] for r1i in r1is_assigned],
        "ais": [[str(coef) for coef in ai_in_p.coefficients] for ai_in_p in ais_in_p],
        "ct0is": [[str(coef) for coef in ct0i_in_p.coefficients] for ct0i_in_p in ct0is_in_p],
    }

    # Calculate the bit size of the largest qi in qis for the filename
    qis_bitsize = max(qis).bit_length()
    qis_len = len(qis)

    # Construct the dynamic filename
    filename = f"sk_enc_{args.n}_{qis_len}x{qis_bitsize}_{args.t}.json"

    output_path = os.path.join("src", "data", "sk_enc_data",filename)

    with open(output_path, 'w') as f:
        json.dump(json_input, f)

    # Initialize a structure to hold polynomials with zero coefficients. This will be used at key generation.
    json_input_zeroes = {
        "s": ["0" for _ in s_assigned.coefficients],
        "e": ["0" for _ in e_assigned.coefficients],
        "k1": ["0" for _ in k1_assigned.coefficients],
        "r2is": [["0" for _ in r2i.coefficients] for r2i in r2is_assigned],
        "r1is": [["0" for _ in r1i.coefficients] for r1i in r1is_assigned],
        "ais": [["0" for _ in ai_in_p.coefficients] for ai_in_p in ais_in_p],
        "ct0is": [["0" for _ in ct0i_in_p.coefficients] for ct0i_in_p in ct0is_in_p],
    }

    output_path = os.path.join("src", "data","sk_enc_data", f"sk_enc_{args.n}_{qis_len}x{qis_bitsize}_{args.t}_zeroes.json")

    with open(output_path, 'w') as f:
        json.dump(json_input_zeroes, f)

    output_path = os.path.join("src", "constants","sk_enc_constants", f"sk_enc_constants_{args.n}_{qis_len}x{qis_bitsize}_{args.t}.rs")

    with open(output_path, 'w') as f:
        f.write(f"/// `N` is the degree of the cyclotomic polynomial defining the ring `Rq = Zq[X]/(X^N + 1)`.\n")
        f.write(f"pub const N: usize = {n};\n")
        f.write(f"/// The coefficients of the polynomial `e` should exist in the interval `[-E_BOUND, E_BOUND]` where `E_BOUND` is the upper bound of the gaussian distribution with 𝜎 = 3.2\n")
        f.write(f"pub const E_BOUND: u64 = {b};\n")
        f.write(f"/// The coefficients of the polynomial `s` should exist in the interval `[-S_BOUND, S_BOUND]`.\n")
        f.write(f"pub const S_BOUND: u64 = {1};\n")
        f.write(f"/// The coefficients of the polynomials `r1is` should exist in the interval `[-R1_BOUND[i], R1_BOUND[i]]` where `R1_BOUND[i]` is equal to `(qi-1)/2`\n")
        f.write(f"pub const R1_BOUNDS: [u64; {len(r1_bounds)}] = [{', '.join(map(str, r1_bounds))}];\n")
        f.write(f"/// The coefficients of the polynomials `r2is` should exist in the interval `[-R2_BOUND[i], R2_BOUND[i]]` where `R2_BOUND[i]` is equal to $\\frac{{(N+2) \\cdot \\frac{{q_i - 1}}{{2}} + B + \\frac{{t - 1}}{{2}} \\cdot |K_{{0,i}}|}}{{q_i}}$\n")
        f.write(f"pub const R2_BOUNDS: [u64; {len(r2_bounds)}] = [{', '.join(map(str, r2_bounds))}];\n")
        f.write(f"/// The coefficients of `k1` should exist in the interval `[-K1_BOUND, K1_BOUND]` where `K1_BOUND` is equal to `(t-1)/2`\n")
        f.write(f"pub const K1_BOUND: u64 = {k1_bound};\n")
        qis_str = ', '.join(f'"{q}"' for q in qi_constants)
        f.write(f"/// List of scalars `qis` such that `qis[i]` is the modulus of the i-th CRT basis of `q` (ciphertext space modulus)\n")
        f.write(f"pub const QIS: [&str; {len(qi_constants)}] = [{qis_str}];\n")
        k0is_str = ', '.join(f'"{k0i}"' for k0i in k0i_constants)
        f.write(f"/// List of scalars `k0is` such that `k0i[i]` is equal to the negative of the multiplicative inverses of t mod qi.\n")
        f.write(f"pub const K0IS: [&str; {len(k0i_constants)}] = [{k0is_str}];\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate rust constants and json inputs for BFV zk proof of secret key encryption circuit"
    )
    parser.add_argument(
        "-n", type=int, required=True, help="Degree of f(x), must be a power of 2."
    )
    parser.add_argument(
        "-qis", type=str, required=True, help="List of qis such that qis[i] is the modulus of the i-th CRT basis of the modulus q of the ciphertext space."
    )
    parser.add_argument(
        "-t", type=int, required=True, help="Modulus t of the plaintext space."
    )

    args = parser.parse_args()
    main(args)
