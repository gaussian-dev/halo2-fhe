import os
from bfv.crt import CRTModuli
from bfv.bfv import BFVCrt
import numpy as np
from bfv.discrete_gauss import DiscreteGaussian
from bfv.polynomial import Polynomial, poly_div
from random import randint
import copy
from utils import assign_to_circuit, count_advice_cells_needed_for_poly_range_check, print_advice_cells_info
import argparse
import json

def main():
    qis =  [1152921504606584833,
            1152921504598720513,
            1152921504597016577,
            1152921504595968001,
            1152921504595640321,
            1152921504593412097,
            1152921504592822273,
            1152921504592429057,
            1152921504589938689,
            1152921504586530817,
            1152921504585547777,
            1152921504583647233,
            1152921504581877761,
            1152921504581419009,
            1152921504580894721]
    crt_moduli = CRTModuli(qis)
    n = 1024
    sigma = 3.2
    discrete_gaussian = DiscreteGaussian(sigma)
    t = 65537
    bfv_crt = BFVCrt(crt_moduli, n, t, discrete_gaussian)

    s = bfv_crt.SecretKeyGen()
    print(f"SECRET KEY LENGTH {len(s.coefficients)}")
    e = bfv_crt.bfv_q.rlwe.SampleFromErrorDistribution()
    ais = []
    for i in range(len(crt_moduli.qis)):
        ais.append(bfv_crt.bfv_qis[i].rlwe.Rq.sample_polynomial())

    pub_keys = bfv_crt.PublicKeyGen(s,e,ais)
    message = bfv_crt.bfv_q.rlwe.Rt.sample_polynomial()

    pk_array = []
    print(f"NUMBER OF PUBLIC KEY {len(pub_keys)}")
   
    for i,pk in enumerate(pub_keys):
        print(f"pk length ={len(pk[0].coefficients)}")
        pk_array.append(pk[0].coefficients)

    print(f"LENGTH OF PK ARRAY {len(pk_array)}")
    


    e0 = bfv_crt.bfv_q.rlwe.SampleFromErrorDistribution()
    e1 = bfv_crt.bfv_q.rlwe.SampleFromErrorDistribution()
    u = bfv_crt.bfv_q.rlwe.SampleFromTernaryDistribution()

    ciphertexts = bfv_crt.PubKeyEncrypt(pub_keys,message,e0,e1,u)
    print(f"NUMBER OF CIPHER TEXTS {len(ciphertexts)}")
    message_prime = bfv_crt.DecryptDummy(s,ciphertexts)




    message_prime == message

    cyclo = [1] + [0] * (n - 1) + [1]
    cyclo = Polynomial(cyclo)

    k1 = message.scalar_mul(crt_moduli.q)
    k1.reduce_coefficients_by_modulus(t)

    for i, cti in enumerate(ciphertexts):
        ct0i =cti[0]

        k0i = pow(-t,-1,qis[i])
        pk0 =  Polynomial(pk_array[i])

        ct0i_hat = pk0 * u + e0 + k1.scalar_mul(k0i)
        assert(len(ct0i_hat.coefficients) - 1 == 2 * n - 2)

        ct0i_hat_clone = copy.deepcopy(ct0i_hat)

        ct0i_hat_clone.reduce_coefficients_by_cyclo(cyclo.coefficients)
        ct0i_hat_clone.reduce_coefficients_by_modulus(qis[i])
        print(f"for iteration ith ={i}")
        assert ct0i_hat_clone == ct0i

    print(" CHECK PASSED ")





main()
    