import os
import argparse
import json
from random import randint
from bfv.crt import CRTModuli
from bfv.bfv import BFVCrt
from bfv.discrete_gauss import DiscreteGaussian
from bfv.polynomial import Polynomial
from utils import assign_to_circuit, count_advice_cells_needed_for_poly_range_check, print_advice_cells_info,poly_range_check,CiphertextFormationInput,range_check_ciphertext_formation
from rlwe import rlwe_encryption


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

    # Generate Public key
    pub_key = bfv_crt.PublicKeyGen(s,e,ais)

    m = bfv_crt.bfv_q.rlwe.Rt.sample_polynomial()
    e0 = bfv_crt.bfv_q.rlwe.SampleFromErrorDistribution()
    e1 = bfv_crt.bfv_q.rlwe.SampleFromErrorDistribution()
    u = bfv_crt.bfv_q.rlwe.SampleFromTernaryDistribution()


    ciphertext = bfv_crt.PubKeyEncrypt(pub_key,m,e0,e1,u)

    # Sanity check for valid decryption 
    message_prime = bfv_crt.Decrypt(s, ciphertext)

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

    # `ct1is` are the polynomials ct0i for each CRT basis.
    ct1is = []

    # `ct1is_hat` are the polynomials ct0i_hat for each CRT basis.
    ct1is_hat = []

    # `pk0is` are the polynomials pk0 for each i-th CRT basis
    pk0is = []

    # `pk1is` are the polynomials pk0 for each i-th CRT basis
    pk1is = []

    # `p1is` are the polynomials p1i for each i-th CRT basis.
    p1is = []

    # `p2is` are the polynomials p2i for each i-th CRT basis.
    p2is = []


    '''
    SETUP PHASE - performed outside the circuit
    For each CRT basis, we need to compute the polynomials r1i and r2i (check this doc for more details: https://hackmd.io/@gaussian/r1W98Kqqa)
    '''
    pk0_array = []
    pk1_array = []


    for i,pk in enumerate(pub_key):
        pk0_array.append(pk[0])
        pk1_array.append(pk[1])

    cyclo = [1] + [0] * (n - 1) + [1]
    cyclo = Polynomial(cyclo)

    zero_poly = [0] * (n + 1)
    zero_poly = Polynomial(zero_poly)

    for i, cti in enumerate(ciphertext):

        ct0i = cti[0]
        ct1i = cti[1]

        pk0 = pk0_array[i]
        pk1 = pk1_array[i]

        #generates the parameters which are used to calculate ct0i
        ct0i_elements = rlwe_encryption(ct0i, qis[i], u, pk0, e0, k1, cyclo, t, n)

        #generates the parameters which are used to calculate ct1i
        ct1i_elements = rlwe_encryption(ct1i, qis[i], u, pk1, e1, zero_poly, cyclo, t, n)

        pk1is.append(pk1)
        p2is.append(ct1i_elements[0])
        p1is.append(ct1i_elements[1])
        ct1is.append(ct1i_elements[3])
        ct1is_hat.append(ct1i_elements[4])

        pk0is.append(pk0)
        r2is.append(ct0i_elements[0])
        r1is.append(ct0i_elements[1])
        k0is.append(ct0i_elements[2])
        ct0is.append(ct0i_elements[3])
        ct0is_hat.append(ct0i_elements[4])
        
    # `r1_bounds` are the bounds for the coefficients of r1i for each CRT basis
    r1_bounds = []

    # `r2_bounds` are the bounds for the coefficients of r2i for each CRT basis
    r2_bounds = []

    # `p1_bounds` are the bounds for the coefficients of p1i for each CRT basis
    p1_bounds = []

    # `p2_bounds` are the bounds for the coefficients of p2i for each CRT basis
    p2_bounds = []

    # initiate counters for the number of advice cells needed for each constraint phase
    phase_0_assignment_advice_cell_count = 0
    phase_1_range_check_advice_cell_count = 0
    phase_1_eval_at_gamma_constraint_advice_cell_count = 0
    phase_1_encryption_constraint_advice_cell_count = 0

    '''
    CIRCUIT - PHASE 0 - ASSIGNMENT
    '''

    # Every assigned value must be an element of the field Zp. Negative coefficients `-z` are assigned as `p - z`
    pk0i_assigned = []
    pk1i_assigned = []

    for i,pk0i in enumerate(pk0is):
        pk0i_assigned.append(assign_to_circuit(pk0i,p))
        pk1i_assigned.append(assign_to_circuit(pk1is[i],p))
        phase_0_assignment_advice_cell_count += len(pk1i_assigned[i].coefficients)
        phase_0_assignment_advice_cell_count += len(pk0i_assigned[i].coefficients)
    

    e0_assigned = assign_to_circuit(e0, p)
    e1_assigned = assign_to_circuit(e1,p)
    k1_assigned = assign_to_circuit(k1, p)
    u_assigned = assign_to_circuit(u,p)

    phase_0_assignment_advice_cell_count += len(e0_assigned.coefficients)
    phase_0_assignment_advice_cell_count += len(e1_assigned.coefficients)
    phase_0_assignment_advice_cell_count += len(u_assigned.coefficients)
    phase_0_assignment_advice_cell_count += len(k1_assigned.coefficients)

    r1is_assigned = []
    r2is_assigned = []
    p1is_assigned = []
    p2is_assigned = []
    ct0is_assigned = []
    ct1is_assigned = []

    for i in range(len(ciphertext)):
        r1i_assigned = assign_to_circuit(r1is[i], p)
        r2i_assigned = assign_to_circuit(r2is[i], p)
        p1i_assigned = assign_to_circuit(p1is[i],p)
        p2i_assigned = assign_to_circuit(p2is[i],p)
        ct0_assigned = assign_to_circuit(ct0is[i],p)
        ct1_assigned = assign_to_circuit(ct1is[i],p)
        p1is_assigned.append(p1i_assigned)
        p2is_assigned.append(p2i_assigned)
        r1is_assigned.append(r1i_assigned)
        r2is_assigned.append(r2i_assigned)
        ct0is_assigned.append(ct0_assigned)
        ct1is_assigned.append(ct1_assigned)

        phase_0_assignment_advice_cell_count += len(r1i_assigned.coefficients)
        phase_0_assignment_advice_cell_count += len(r2i_assigned.coefficients)
        phase_0_assignment_advice_cell_count += len(p1i_assigned.coefficients)
        phase_0_assignment_advice_cell_count += len(p2i_assigned.coefficients)
        phase_0_assignment_advice_cell_count += len(ct0_assigned.coefficients)
        phase_0_assignment_advice_cell_count += len(ct1_assigned.coefficients)



    # For the sake of simplicity, we generate a random challenge here
    gamma = randint(0, 1000)

    '''
    CIRCUIT - PHASE 1 - ASSIGNMENT
    '''

    # Every assigned value must be an element of the field Zp. Negative coefficients `-z` are assigned as `p - z`
    qi_constants = []
    k0i_constants = []



    for i in range(len(ciphertext)):

        qi_constants.append(qis[i])

        k0i_constant = assign_to_circuit(Polynomial([k0is[i]]), p).coefficients[0]
        k0i_constants.append(k0i_constant)

    cyclo_at_gamma = cyclo.evaluate(gamma)
    cyclo_at_gamma_assigned = assign_to_circuit(Polynomial([cyclo_at_gamma]), p).coefficients[0]


    '''
    CIRCUIT - PHASE 1 - RANGE CHECK
    '''

    lookup_bits = 8

    b = int(discrete_gaussian.z_upper)
    # constraint. The coefficients of e should be in the range [-B, B] where B is the upper bound of the discrete Gaussian distribution
    poly_range_check(e0,e0_assigned,p,b)

    phase_1_range_check_advice_cell_count += count_advice_cells_needed_for_poly_range_check(e0_assigned, 2*b + 1, lookup_bits)

    # constraint. The coefficients of e1 should be in the range [-B, B] where B is the upper bound of the discrete Gaussian distribution
    poly_range_check(e1,e1_assigned,p,b)
    phase_1_range_check_advice_cell_count += count_advice_cells_needed_for_poly_range_check(e1_assigned, 2*b + 1, lookup_bits)

    u_bound = 1
    # constraint. The coefficients of u should be in the range [-1, 0, 1]
    poly_range_check(u,u_assigned,p,u_bound)
    phase_1_range_check_advice_cell_count += count_advice_cells_needed_for_poly_range_check(u_assigned, 2*u_bound + 1, lookup_bits)

    k1_bound = int((t - 1) / 2)
    # constraint. The coefficients of k1 should be in the range [-(t-1)/2, (t-1)/2]
    poly_range_check(k1,k1_assigned,p,k1_bound)

    phase_1_range_check_advice_cell_count += count_advice_cells_needed_for_poly_range_check(k1_assigned, 2*k1_bound + 1, lookup_bits)


    e0_at_gamma_assigned = e0_assigned.evaluate(gamma)
    phase_1_eval_at_gamma_constraint_advice_cell_count += len(e0_assigned.coefficients) * 2 - 1

    e1_at_gamma_assigned = e1_assigned.evaluate(gamma)
    phase_1_eval_at_gamma_constraint_advice_cell_count += len(e1_assigned.coefficients) * 2 - 1


    k1_at_gamma_assigned = k1_assigned.evaluate(gamma)
    phase_1_eval_at_gamma_constraint_advice_cell_count += len(k1_assigned.coefficients) * 2 - 1

    u_at_gamma_assigned = u_assigned.evaluate(gamma)
    phase_1_eval_at_gamma_constraint_advice_cell_count += len(u_assigned.coefficients) * 2 - 1


    pk_bound = []

    for i in range(len(ciphertext)):
        args = CiphertextFormationInput(
            qis=qis[i],
            n=n,
            b=b,
            t=t,
            ctis=ct0is[i],
            ais_pk0=pk0is[i],
            s_u=u,
            e=e0,
            r2is=r2is[i],
            k0is=k0is[i],
            k1=k1,
            ctis_hat=ct0is_hat[i],
            r1is=r1is[i]
        )

        #conducts a range check at each step of forming ct0i
        range_check_ciphertext_formation(args)

        #constraint .The coefficient of pk0i_assigned and pk1i_assigned should be in range [-(qi-1)/2 , (qi-1)/2 ]
        pk0_bound = int((qis[i] - 1) / 2)
        pk_bound.append(pk0_bound)
        poly_range_check(pk0is[i],pk0i_assigned[i],p,pk0_bound)
        poly_range_check(pk1is[i],pk1i_assigned[i],p,pk0_bound)


        phase_1_range_check_advice_cell_count += count_advice_cells_needed_for_poly_range_check(pk0i_assigned[i], 2*pk0_bound + 1, lookup_bits)
        phase_1_range_check_advice_cell_count += count_advice_cells_needed_for_poly_range_check(pk1i_assigned[i], 2*pk0_bound + 1, lookup_bits)


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

        # constraint  the coefficients of p2 should be in the range [-(qi-1)/2, (qi-1)/2]
        # p2i should be in the modulus so the range of p1i should be qi - 1 ,but during range check we add the upperbound to the element to control the negative elements
        # so, the p2i should be in the range (qi - 1) / 2 ,such that 2 * p2i is within qi - 1 
        p2i_bound = int((qis[i] - 1) / 2)
        p2_bounds.append(p2i_bound)
        poly_range_check(p2is[i],p2is_assigned[i],p,p2i_bound)

        phase_1_range_check_advice_cell_count += count_advice_cells_needed_for_poly_range_check(p2is_assigned[i], 2*p2i_bound + 1, lookup_bits)

        #ct0i - ct0i_hat = $[ (N+1) \cdot \frac{q_i - 1}{2} + B , (N+1) \cdot \frac{q_i - 1}{2} + B ]$
        #ct0i - ct0i_hat - p2i * cyclo = $[ (N+2) \cdot \frac{q_i - 1}{2} + B , (N+2) \cdot \frac{q_i - 1}{2} + B ]$
        # p1i = (ct0i - ct0i_hat - p2i * cyclo)/ qi = $[\frac{- ((N+2) \cdot \frac{q_i - 1}{2} + B )}{q_i}, \frac{((N+2) \cdot \frac{q_i - 1}{2} + B )}{q_i}]$
        # constraint. The coefficients of (ct0i - ct0i_hat - p2i * cyclo) / qi = p1i should be in the range [(((qis[i] - 1) / 2) * (n + 2) + b ) /qis[i]]
        p1i_bound = (int((qis[i] - 1) / 2) * (n + 2) + b ) / qis[i]
        # round bound to the nearest integer
        p1i_bound = int(p1i_bound)
        p1_bounds.append(p1i_bound)
        poly_range_check(p1is[i],p1is_assigned[i],p,p1i_bound)

        phase_1_range_check_advice_cell_count += count_advice_cells_needed_for_poly_range_check(p1is_assigned[i], 2*p1i_bound + 1, lookup_bits)

        # conducts a range check at each step of forming ct1i
        # sanity check. The coefficients of p2is[i] * cyclo should be in the range [-(qi-1)/2, (qi-1)/2]
        bound = int((qis[i] - 1) / 2)
        res = p2is[i] * cyclo
        assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

        # sanity check. The coefficients of ct1i[i] - ct1i_hat[i] should be in the range $ [- ((N+1) \cdot \frac{q_i - 1}{2} + B), (N+1) \cdot \frac{q_i - 1}{2} + B]$
        bound = int((qis[i] - 1) / 2) * (n + 1) + b
        sub = ct1is[i] + ct1is_hat[i].scalar_mul(-1)
        assert all(coeff >= -bound and coeff <= bound for coeff in sub.coefficients)

        # sanity check. The coefficients of ct1i[i] - ct1i_hat[i] - p2is[i] * cyclo should be in the range $[- ((N+2) \cdot \frac{q_i - 1}{2} + B), (N+2) \cdot \frac{q_i - 1}{2} + B]$
        bound = int((qis[i] - 1) / 2) * (n + 2) + b
        sub = ct1is[i] + (ct1is_hat[i].scalar_mul(-1)) +  (p2is[i] * cyclo).scalar_mul(-1)
        assert all(coeff >= -bound and coeff <= bound for coeff in sub.coefficients)


        '''
        CIRCUIT - PHASE 1 - EVALUATION AT GAMMA CONSTRAINT
        '''

        r1i_gamma_assigned = r1is_assigned[i].evaluate(gamma)
        r2i_gamma_assigned = r2is_assigned[i].evaluate(gamma)

        phase_1_eval_at_gamma_constraint_advice_cell_count += len(r1is_assigned[i].coefficients) * 2 - 1
        phase_1_eval_at_gamma_constraint_advice_cell_count += len(r2is_assigned[i].coefficients) * 2 - 1

        pk0i_at_gamma = pk0is[i].evaluate(gamma)
        pk0i_at_gamma_assigned = assign_to_circuit(Polynomial([pk0i_at_gamma]),p).coefficients[0]

        p1i_gamma_assigned = p1is_assigned[i].evaluate(gamma)
        p2i_gamma_assigned = p2is_assigned[i].evaluate(gamma)

        phase_1_eval_at_gamma_constraint_advice_cell_count += len(p1is_assigned[i].coefficients) * 2 - 1
        phase_1_eval_at_gamma_constraint_advice_cell_count += len(p2is_assigned[i].coefficients) * 2 - 1

        pk1i_at_gamma = pk1is[i].evaluate(gamma)
        pk1i_at_gamma_assigned = assign_to_circuit(Polynomial([pk1i_at_gamma]),p).coefficients[0]


        '''
        CIRCUIT - PHASE 1 - CORRECT ENCRYPTION CONSTRAINT
        '''

        lhs = ct0is_assigned[i].evaluate(gamma)
        rhs = (pk0i_at_gamma_assigned * u_at_gamma_assigned + e0_at_gamma_assigned + (k1_at_gamma_assigned * k0i_constants[i]) + (r1i_gamma_assigned * qi_constants[i]) + (r2i_gamma_assigned * cyclo_at_gamma_assigned))
        phase_1_encryption_constraint_advice_cell_count += 32

        assert lhs % p == rhs % p

        lhs = ct1is_assigned[i].evaluate(gamma)
        rhs = (pk1i_at_gamma_assigned * u_at_gamma_assigned + e1_at_gamma_assigned + p2i_gamma_assigned * cyclo_at_gamma_assigned + p1i_gamma_assigned * qi_constants[i])

        assert lhs % p == rhs % p
        '''
        VERIFICATION PHASE
        '''

        cyclo_at_gamma = cyclo.evaluate(gamma)
        cyclo_at_gamma_assigned_expected = assign_to_circuit(Polynomial([cyclo_at_gamma]), p).coefficients[0]
        assert cyclo_at_gamma_assigned == cyclo_at_gamma_assigned_expected

        
        pk0i_gamma = pk0is[i].evaluate(gamma)
        pk0i_gamma_assigned_expected = assign_to_circuit(Polynomial([pk0i_gamma]), p).coefficients[0]
        assert pk0i_at_gamma_assigned == pk0i_gamma_assigned_expected


        pk1i_gamma = pk1is[i].evaluate(gamma)
        pk1i_gamma_assigned_expected = assign_to_circuit(Polynomial([pk1i_gamma]), p).coefficients[0]
        assert pk1i_at_gamma_assigned == pk1i_gamma_assigned_expected 



        assert qis[i] == qi_constants[i]

        k0i_assigned_expected = assign_to_circuit(Polynomial([k0is[i]]), p).coefficients[0]
        assert k0i_constants[i] == k0i_assigned_expected

    total_advice_cell_count = phase_0_assignment_advice_cell_count + phase_1_range_check_advice_cell_count + phase_1_eval_at_gamma_constraint_advice_cell_count + phase_1_encryption_constraint_advice_cell_count

    print_advice_cells_info(total_advice_cell_count, phase_0_assignment_advice_cell_count, phase_1_range_check_advice_cell_count, phase_1_eval_at_gamma_constraint_advice_cell_count, phase_1_encryption_constraint_advice_cell_count)
    #  ct0is and ct1is need to be parsed such that their coefficients are in the range [0, p - 1]
    ct0is_in_p = [assign_to_circuit(ct0i, p) for ct0i in ct0is]
    ct1is_in_p = [assign_to_circuit(ct1i, p) for ct1i in ct1is]

    # Parse the inputs into a JSON format such this can be used as input for the (real) circuit
    json_input = {
        "pk0_qi": [[str(coef) for coef in pk0i.coefficients] for pk0i in pk0i_assigned],
        "pk1_qi": [[str(coef) for coef in pk1i.coefficients] for pk1i in pk1i_assigned],
        "u": [str(coef) for coef in u_assigned.coefficients],
        "e0": [str(coef) for coef in e0_assigned.coefficients],
        "e1": [str(coef) for coef in e1_assigned.coefficients],
        "k1": [str(coef) for coef in k1_assigned.coefficients],
        "r2is": [[str(coef) for coef in r2i.coefficients] for r2i in r2is_assigned],
        "r1is": [[str(coef) for coef in r1i.coefficients] for r1i in r1is_assigned],
        "p2is": [[str(coef) for coef in p2i.coefficients] for p2i in p2is_assigned],
        "p1is": [[str(coef) for coef in p1i.coefficients] for p1i in p1is_assigned],
        "ct0is": [[str(coef) for coef in ct0i_in_p.coefficients] for ct0i_in_p in ct0is_in_p],
        "ct1is": [[str(coef) for coef in ct0i_in_p.coefficients] for ct0i_in_p in ct1is_in_p],
    }

    # Calculate the bit size of the largest qi in qis for the filename
    qis_bitsize = max(qis).bit_length()
    qis_len = len(qis)

    # Construct the dynamic filename
    filename = f"pk_enc_{args.n}_{qis_len}x{qis_bitsize}_{args.t}.json"

    output_path = os.path.join("src", "data","pk_enc_data", filename)

    with open(output_path, 'w') as f:
        json.dump(json_input, f)

    # Initialize a structure to hold polynomials with zero coefficients. This will be used at key generation.
    json_input_zeroes = {
        "pk0_qi":[["0" for _ in pk0i.coefficients] for pk0i in pk0i_assigned],
        "pk1_qi":[["0" for _ in pk1i.coefficients] for pk1i in pk1i_assigned],
        "u": ["0" for _ in u_assigned.coefficients],
        "e0": ["0" for _ in e0_assigned.coefficients],
        "e1": ["0" for _ in e1_assigned.coefficients],
        "k1": ["0" for _ in k1_assigned.coefficients],
        "r2is": [["0" for _ in r2i.coefficients] for r2i in r2is_assigned],
        "r1is": [["0" for _ in r1i.coefficients] for r1i in r1is_assigned],
        "p2is": [["0" for _ in p2i.coefficients] for p2i in p2is_assigned],
        "p1is": [["0" for _ in p1i.coefficients] for p1i in p1is_assigned],
        "ct0is": [["0" for _ in ct0i_in_p.coefficients] for ct0i_in_p in ct0is_in_p],
        "ct1is": [["0" for _ in ct1i_in_p.coefficients] for ct1i_in_p in ct1is_in_p],
    }

    output_path = os.path.join("src", "data","pk_enc_data", f"pk_enc_{args.n}_{qis_len}x{qis_bitsize}_{args.t}_zeroes.json")

    with open(output_path, 'w') as f:
        json.dump(json_input_zeroes, f)

    output_path = os.path.join("src", "constants","pk_enc_constants", f"pk_enc_constants_{args.n}_{qis_len}x{qis_bitsize}_{args.t}.rs")

    with open(output_path, 'w') as f:
        f.write(f"/// `N` is the degree of the cyclotomic polynomial defining the ring `Rq = Zq[X]/(X^N + 1)`.\n")
        f.write(f"pub const N: usize = {n};\n")
        f.write(f"///'The coefficients pf the polynomial 'pk0is` and 'pk1is' should exist in the interval '[-PK_BOUND, PK_BOUND]`.\n")
        pk_bound_str = ', '.join(map(str,pk_bound))
        f.write(f"pub const PK_BOUND :[u64; {len(pk_bound)}] = [{pk_bound_str}];\n")
        f.write(f"///'The coefficients pf the polynomial 'pk1is` should exist in the interval '[-PK0_BOUND, PK0_BOUND]`.\n")
        f.write(f"/// The coefficients of the polynomial `e` should exist in the interval `[-E_BOUND, E_BOUND]` where `E_BOUND` is the upper bound of the gaussian distribution with 𝜎 = 3.2\n")
        f.write(f"pub const E_BOUND: u64 = {b};\n")
        f.write(f"/// The coefficients of the polynomial `s` should exist in the interval `[-S_BOUND, S_BOUND]`.\n")
        f.write(f"pub const U_BOUND: u64 = {1};\n")
        f.write(f"/// The coefficients of the polynomials `r1is` should exist in the interval `[-R1_BOUND[i], R1_BOUND[i]]` where `R1_BOUND[i]` is equal to `(qi-1)/2`\n")
        f.write(f"pub const R1_BOUNDS: [u64; {len(r1_bounds)}] = [{', '.join(map(str, r1_bounds))}];\n")
        f.write(f"/// The coefficients of the polynomials `r2is` should exist in the interval `[-R2_BOUND[i], R2_BOUND[i]]` where `R2_BOUND[i]` is equal to $\\frac{{(N+2) \\cdot \\frac{{q_i - 1}}{{2}} + B + \\frac{{t - 1}}{{2}} \\cdot |K_{{0,i}}|}}{{q_i}}$\n")
        f.write(f"pub const R2_BOUNDS: [u64; {len(r2_bounds)}] = [{', '.join(map(str, r2_bounds))}];\n")
        f.write(f"/// The coefficients of the polynomials `p1is` should exist in the interval `[-P1_BOUND[i], P1_BOUND[i]]` where `P1_BOUND[i]` is equal to (((qis[i] - 1) / 2) * (n + 2) + b ) / qis[i] \n")
        f.write(f"pub const P1_BOUNDS: [u64; {len(p1_bounds)}] = [{', '.join(map(str, p1_bounds))}];\n")
        f.write(f"/// The coefficients of the polynomials `p2is` should exist in the interval `[-P2_BOUND[i], P2_BOUND[i]]` where `P2_BOUND[i]` is equal to (qis[i] - 1) / 2  \n")
        f.write(f"pub const P2_BOUNDS: [u64; {len(p2_bounds)}] = [{', '.join(map(str, p2_bounds))}];\n")
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
