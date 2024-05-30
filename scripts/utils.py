import math
import copy
from bfv.polynomial import Polynomial,poly_div


def rlwe_encryption(cti: Polynomial, qi: Polynomial, u_s: Polynomial, pk_a: Polynomial, e: Polynomial, k1: Polynomial, cyclo: Polynomial, t: int, n: int) -> list:
    '''
    This function takes input and performs the rlwe encryption 
    `cti` is ciphertext of either pk-enc or sk-enc
    `u_s` it contains either `u` or `s` based on the encryption used
    `pk_a` it contains either `pk0`, `pk1` or `a` based on the encryption used
    '''
    enc_elements = []

    # k0i = -t^{-1} namely the multiplicative inverse of t modulo qi
    k0i = pow(-t,-1,qi)

    cti_hat = pk_a * u_s + e + k1.scalar_mul(k0i)
    assert(len(cti_hat.coefficients) - 1 == 2 * n - 2)
    # pk_a * u_s + e + k0i * k1 = cti mod Rqi
    # assert that cti_hat = cti mod Rqi
    cti_hat_clone = copy.deepcopy(cti_hat)
    # mod Rqi means that we need to:
    # - reduce the coefficients of cti_hat_clone by the cyclotomic polynomial
    # - reduce the coefficients of cti_hat_clone by the modulus
    cti_hat_clone.reduce_coefficients_by_cyclo(cyclo.coefficients)
    cti_hat_clone.reduce_coefficients_by_modulus(qi)
    assert cti_hat_clone == cti

    
    # Calculate r2i_p2i 
    # divide cti - cti_hat by the cyclotomic polynomial over Zqi to get r2i or p2i based on the encryption used.
    num = cti + cti_hat.scalar_mul(-1)
    # reduce the coefficients of num by the modulus qi 
    num.reduce_coefficients_by_modulus(qi)
    (quotient, rem) = poly_div(num.coefficients, cyclo.coefficients)
    # assert that the remainder is zero
    assert rem == []
    r2i_p2i = Polynomial(quotient)
    # assert that the degree of r2i_p2i is equal to n - 2
    assert len(r2i_p2i.coefficients) - 1 == n - 2
    enc_elements.append(r2i_p2i)

    # Assert that cti - cti_hat = r2i_p2i * cyclo mod Zqi
    lhs = cti + cti_hat.scalar_mul(-1)
    rhs = r2i_p2i * cyclo
    # reduce the coefficients of lhs by the modulus qi
    lhs.reduce_coefficients_by_modulus(qi)
    assert lhs == rhs
        
    # Calculate r1i
    # divide cti - cti_hat - r2i_p2i * cyclo by the modulus qi to get r1i
    num = cti + cti_hat.scalar_mul(-1) + (r2i_p2i * cyclo).scalar_mul(-1)
    (quotient, rem) = poly_div(num.coefficients, [qi])
    # assert that the remainder is zero
    assert rem == []
    r1i = Polynomial(quotient)
    # assert that the degree of r1i is 2n - 2
    assert len(r1i.coefficients) - 1 == 2 * n - 2
    enc_elements.append(r1i)

    enc_elements.append(k0i)

    # Assert that cti = cti_hat + r1i * qi + r2i_p2i * cyclo mod Zp
    lhs = cti
    enc_elements.append(cti)

    rhs = cti_hat + (r1i.scalar_mul(qi)) + (r2i_p2i * cyclo)
    enc_elements.append(cti_hat)

    # remove the leading zeroes from rhs until the length of rhs.coefficients is equal to n
    while len(rhs.coefficients) > n and rhs.coefficients[0] == 0:
        rhs.coefficients.pop(0)


    assert lhs == rhs
    return enc_elements


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
