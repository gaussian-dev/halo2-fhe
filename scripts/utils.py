import math
from typing import NamedTuple
from bfv.polynomial import Polynomial

class CiphertextFormationInput(NamedTuple):
    '''
    This class includes all the inputs required for the calculation of the ciphertext(1st Matrix)
    'ais_pk0' is a polynomial, with 'ais' used for sk_enc and `pk0` used for pk_enc
    `s_u` is a polynomial, with 'ais' used for sk_enc and `pk0` used for pk_enc
    '''
    qis: int
    n: int
    b: int
    t: int
    ctis: Polynomial
    ais_pk0: Polynomial
    s_u: Polynomial
    e: Polynomial
    r2is: Polynomial
    k0is: Polynomial
    k1:Polynomial
    ctis_hat: Polynomial
    r1is: Polynomial


def poly_range_check(poly:Polynomial , poly_assigned:Polynomial, p:int, bound:int):
    '''
    This function takes `bound` and `polynomial` and checks its shifted coefficients are within 2*bound
    `poly` is the polynomial in its original form
    `poly_assigned` is the polynomial that is assignedd to the circuit 
    `bound` is the range in which the polynimial coefficient should exist 
    '''

    # constraint. The coefficients of poly should be in the range [-bound, 0, bound]
    assert all(coeff >= -bound and coeff <= bound for coeff in poly.coefficients)
    # After the circuit assignement, the coefficients of poly_assigned must be in  [0, bound] or [p - bound, p - 1]
    assert all(coeff in range(0, bound+1) or coeff in range(p - bound, p) for coeff in poly_assigned.coefficients)
    # To perform a range check with a smaller lookup table, we shift the coefficients of poly_assigned to be in [0, 2B] (the shift operation is constrained inside the circuit)
    poly_shifted = Polynomial([(coeff + bound) % p for coeff in poly_assigned.coefficients])
    assert all(coeff >= 0 and coeff <= 2*bound for coeff in poly_shifted.coefficients)


def range_check_ciphertext_formation(args:CiphertextFormationInput):

    '''
    This function takes arguments as input, performs range checks on the polynomial after each operation,
    and calculates the first matrix (ct0) for sk_enc and pk_enc.
    '''

    cyclo = [1] + [0] * (args.n - 1) + [1]
    cyclo = Polynomial(cyclo)
    # sanity check. The coefficients of ct0i should be in the range [-(qi-1)/2, (qi-1)/2]
    bound = int((args.qis - 1) / 2)
    assert all(coeff >= -bound and coeff <= bound for coeff in args.ctis.coefficients)

    # sanity check. The coefficients of `ais` or `pk0` should be in the range [-(qi-1)/2, (qi-1)/2]
    bound = int((args.qis - 1) / 2)
    assert all(coeff >= -bound and coeff <= bound for coeff in args.ais_pk0.coefficients)

    # sanity check. The coefficients of ais_pk0 * s_u should be in the range $[-N \cdot \frac{q_i - 1}{2}, N \cdot \frac{q_i - 1}{2}]$
    bound = int((args.qis - 1) / 2) * args.n
    res = args.ais_pk0 * args.s_u
    assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

    # sanity check. The coefficients of ais_pk0 * s_u + e should be in the range $- (N \cdot \frac{q_i - 1}{2} + B), N \cdot \frac{q_i - 1}{2} + B]$
    bound = int((args.qis - 1) / 2) * args.n + args.b
    res = args.ais_pk0 * args.s_u + args.e
    assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

    # sanity check. The coefficients of r2i * cyclo should be in the range [-(qi-1)/2, (qi-1)/2]
    bound = int((args.qis - 1) / 2)
    res = args.r2is * cyclo
    assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

    # sanity check. The coefficients of k1 * k0i should be in the range $[-\frac{t - 1}{2} \cdot |K_{0,i}|, \frac{t - 1}{2} \cdot |K_{0,i}|]$
    bound = int((args.t - 1) / 2) * abs(args.k0is)
    res = args.k1.scalar_mul(args.k0is)
    assert all(coeff >= -bound and coeff <= bound for coeff in res.coefficients)

    # sanity check. The coefficients of ct0i_hat (ais_pk0 * s_u + e + k1 * k0i) should be in the range $[- (N \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_{0,i}|), N \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_{0,i}|]$
    bound = int((args.qis - 1) / 2) * args.n + args.b + int((args.t - 1) / 2) * abs(args.k0is)
    assert all(coeff >= -bound and coeff <= bound for coeff in args.ctis_hat.coefficients)

    # sanity check. The coefficients of ct0i - ct0i_hat should be in the range $ [- ((N+1) \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_{0,i}|), (N+1) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_{0,i}|]$
    bound = int((args.qis - 1) / 2) * (args.n + 1) + args.b + int((args.t - 1) / 2) * abs(args.k0is)
    sub = args.ctis + (args.ctis_hat.scalar_mul(-1))
    assert all(coeff >= -bound and coeff <= bound for coeff in sub.coefficients)

    # sanity check. The coefficients of ct0i - ct0i_hat - r2i * cyclo should be in the range $[- ((N+2) \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_{0,i}|), (N+2) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_{0,i}|]$
    bound = ((args.qis - 1) / 2) * (args.n + 2) + args.b + ((args.t - 1) / 2) * abs(args.k0is)
    sub = args.ctis + (args.ctis_hat.scalar_mul(-1)) + (args.r2is * cyclo).scalar_mul(-1)
    assert all(coeff >= -bound and coeff <= bound for coeff in sub.coefficients)



def assign_to_circuit(poly: Polynomial, p: int) -> Polynomial:
    '''
    This function takes a polynomial and returns its coefficients in the field Zp
    `poly` is the polynomial to be assigned to the circuit
    `p` is the field modulus
    '''
    assigned_coefficients = []
    for coeff in poly.coefficients:
        if coeff < 0:
            coeff = coeff % p
        if coeff > p:
            coeff = coeff % p
        assigned_coefficients.append(coeff)

    return Polynomial(assigned_coefficients)

def count_advice_cells_needed_for_poly_range_check(poly: Polynomial, bound: int, lookup_bits: int) -> int:
    '''
    This function takes a polynomial and a bound and returns the number of advice cells needed for a range check
    `poly` is the polynomial to be checked
    `bound` is the upper bound for the range check
    `lookup_bits` is the number of bits used for the lookup table
    '''

    count = 0

    # 4 advice cells for each coefficient needed for the shift addition operation``
    count += 4 * len(poly.coefficients)

    # further advice cells for range check inside `check_less_than_safe`
    bound_bits = bound.bit_length()
    range_bits = math.ceil(bound_bits / lookup_bits) * lookup_bits
    num_limbs = math.ceil(range_bits / lookup_bits)

    if num_limbs > 1:
        # 1 + (3 * (num_limbs - 1)) advice cells
        count += (1 + (3 * (num_limbs - 1))) * len(poly.coefficients)
    else:
        # count is not updated if num_limbs is 1
        pass

    # 7 advice cells required for the `check_less_than` constraint inside `check_less_than_safe`
    count += 7 * len(poly.coefficients)

    # the check_less_than_advice_cells constraint also performs a range check on the check_cell in range_bits
    if num_limbs > 1:
        # 1 + (3 * (num_limbs - 1)) advice cells for the check_less_than_advice_cells constraint
        count += (1 + (3 * (num_limbs - 1))) * len(poly.coefficients)
    else:
        # count is not updated if num_limbs is 1
        pass

    return count

def print_advice_cells_info(total_advice_cell_count, phase_0_count, phase_1_range_check_count, phase_1_eval_at_gamma_count, phase_1_encryption_constraint_count):
    print("Halo2 Circuit Profile:")
    print(f"Total Advice Cells Needed: {total_advice_cell_count}")
    
    print("\nPhase 0 - Assignment:")
    print(f" - Count: {phase_0_count}, Percentage: {(phase_0_count / total_advice_cell_count) * 100:.2f}%")

    print("\nPhase 1 - Range Check:")
    print(f" - Count: {phase_1_range_check_count}, Percentage: {(phase_1_range_check_count / total_advice_cell_count) * 100:.2f}%")

    print("\nPhase 1 - Evaluation at Gamma Constraint:")
    print(f" - Count: {phase_1_eval_at_gamma_count}, Percentage: {(phase_1_eval_at_gamma_count / total_advice_cell_count) * 100:.2f}%")

    print("\nPhase 1 - Correct Encryption Constraint:")
    print(f" - Count: {phase_1_encryption_constraint_count}, Percentage: {(phase_1_encryption_constraint_count / total_advice_cell_count) * 100:.2f}%")
