use crate::amcl_wrapper::group_elem::{GroupElement, GroupElementVector};
use amcl_wrapper::extension_field_gt::GT;
use amcl_wrapper::group_elem_g1::{G1, G1Vector};
use amcl_wrapper::group_elem_g2::{G2, G2Vector};
use crate::transcript::TranscriptProtocol;
use merlin::Transcript;
use amcl_wrapper::field_elem::FieldElement;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct InnerPairingProductProof {
    pub T_L: Vec<GT>,
    pub T_R: Vec<GT>,
    pub U_L: Vec<GT>,
    pub U_R: Vec<GT>,
    pub Z_L: Vec<GT>,
    pub Z_R: Vec<GT>,
    pub A: G1,
    pub B: G2
}

pub struct InnerPairingProduct {}

impl InnerPairingProduct {
    pub fn gen_proof(transcript: &mut Transcript, T: &GT, U: &GT, Z: &GT, w: &G1Vector, v: &G2Vector, A: &G1Vector, B: &G2Vector) -> InnerPairingProductProof {
        let mut m = w.len();
        let lg_m = m.next_power_of_two().trailing_zeros() as usize;

        // m is a power of 2. Not doing n & n-1 == 0 check as log_2 m is needed
        assert_eq!(1 << lg_m, m);

        assert_eq!(v.len(), m);
        assert_eq!(A.len(), m);
        assert_eq!(B.len(), m);

        let mut w = w.clone();
        let mut v = v.clone();
        let mut A = A.clone();
        let mut B = B.clone();

        let mut T_L_vec = Vec::<GT>::with_capacity(lg_m);
        let mut T_R_vec = Vec::<GT>::with_capacity(lg_m);
        let mut U_L_vec = Vec::<GT>::with_capacity(lg_m);
        let mut U_R_vec = Vec::<GT>::with_capacity(lg_m);
        let mut Z_L_vec = Vec::<GT>::with_capacity(lg_m);
        let mut Z_R_vec = Vec::<GT>::with_capacity(lg_m);

        while m != 1 {
            m >>= 1;

            let (mut w_L, w_R) = w.split_at(m);
            let (mut v_L, v_R) = v.split_at(m);
            let (mut A_L, A_R) = A.split_at(m);
            let (mut B_L, B_R) = B.split_at(m);

            let T_L = GT::inner_product(A_L.as_slice(), v_R.as_slice()).unwrap();
            let T_R = GT::inner_product(A_R.as_slice(), v_L.as_slice()).unwrap();
            let U_L = GT::inner_product(w_L.as_slice(), B_R.as_slice()).unwrap();
            let U_R = GT::inner_product(w_R.as_slice(), B_L.as_slice()).unwrap();
            let Z_L = GT::inner_product(A_L.as_slice(), B_R.as_slice()).unwrap();
            let Z_R = GT::inner_product(A_R.as_slice(), B_L.as_slice()).unwrap();

            transcript.commit_GT(b"T_L", &T_L);
            transcript.commit_GT(b"T_R", &T_R);
            transcript.commit_GT(b"U_L", &U_L);
            transcript.commit_GT(b"U_R", &U_R);
            transcript.commit_GT(b"Z_L", &Z_L);
            transcript.commit_GT(b"Z_R", &Z_R);

            T_L_vec.push(T_L);
            T_R_vec.push(T_R);
            U_L_vec.push(U_L);
            U_R_vec.push(U_R);
            Z_L_vec.push(Z_L);
            Z_R_vec.push(Z_R);

            let x = transcript.challenge_scalar(b"x");
            let x_inv = x.inverse();

            for i in 0..m {
                // A_L[i] = (x * A_L[i]) + (x_inv * A_R[i])
                A_L[i] = A_L[i].binary_scalar_mul(&A_R[i], &x, &x_inv);
                // B_L[i] = (x_inv * B_L[i]) + (x * B_R[i])
                B_L[i] = B_L[i].binary_scalar_mul(&B_R[i], &x_inv, &x);
                // w_L[i] = (x * w_L[i]) + (x_inv * w_R[i])
                w_L[i] = w_L[i].binary_scalar_mul(&w_R[i], &x, &x_inv);
                // v_L[i] = (x * v_L[i]) + (x_inv * v_R[i])
                v_L[i] = v_L[i].binary_scalar_mul(&v_R[i], &x_inv, &x);
            }

            A = A_L;
            B = B_L;
            w = w_L;
            v = v_L;
        }

        assert_eq!(A.len(), 1);
        assert_eq!(B.len(), 1);
        InnerPairingProductProof {
            T_L: T_L_vec,
            T_R: T_R_vec,
            U_L: U_L_vec,
            U_R: U_R_vec,
            Z_L: Z_L_vec,
            Z_R: Z_R_vec,
            A: A.pop().unwrap(),
            B: B.pop().unwrap(),
        }
    }

    pub fn verify_proof_recursively(m: usize, transcript: &mut Transcript, T: &GT, U: &GT, Z: &GT, w: &G1Vector, v: &G2Vector, proof: &InnerPairingProductProof) -> bool {
        let lg_m = m.next_power_of_two().trailing_zeros() as usize;
        // m is a power of 2. Not doing n & n-1 == 0 check as log_2 m is needed
        assert_eq!(1 << lg_m, m);

        assert_eq!(w.len(), m);
        assert_eq!(v.len(), m);

        assert_eq!(proof.T_L.len(), proof.T_R.len());
        assert_eq!(proof.U_L.len(), proof.U_R.len());
        assert_eq!(proof.Z_L.len(), proof.Z_R.len());
        assert_eq!(proof.T_L.len(), lg_m);
        assert_eq!(proof.U_L.len(), lg_m);
        assert_eq!(proof.Z_L.len(), lg_m);

        let mut m = m;
        let mut T = T.clone();
        let mut U = U.clone();
        let mut Z = Z.clone();
        let mut w = w.clone();
        let mut v = v.clone();

        let mut i = 0;
        while m != 1 {
            m >>= 1;

            let (mut w_L, w_R) = w.split_at(m);
            let (mut v_L, v_R) = v.split_at(m);

            let T_L = &proof.T_L[i];
            let T_R = &proof.T_R[i];
            let U_L = &proof.U_L[i];
            let U_R = &proof.U_R[i];
            let Z_L = &proof.Z_L[i];
            let Z_R = &proof.Z_R[i];

            transcript.commit_GT(b"T_L", &T_L);
            transcript.commit_GT(b"T_R", &T_R);
            transcript.commit_GT(b"U_L", &U_L);
            transcript.commit_GT(b"U_R", &U_R);
            transcript.commit_GT(b"Z_L", &Z_L);
            transcript.commit_GT(b"Z_R", &Z_R);

            let x = transcript.challenge_scalar(b"x");
            let x_inv = x.inverse();
            let x_sqr = x.square();
            // Square is cheaper than inverse so not doing x_sqr.inverse()
            let x_inv_sqr = x_inv.square();

            for i in 0..m {
                // w_L[i] = (x * w_L[i]) + (x_inv * w_R[i])
                w_L[i] = w_L[i].binary_scalar_mul(&w_R[i], &x, &x_inv);
                // v_L[i] = (x * v_L[i]) + (x_inv * v_R[i])
                v_L[i] = v_L[i].binary_scalar_mul(&v_R[i], &x_inv, &x);
            }

            T = (T_L.pow(&x_sqr) * &T) * T_R.pow(&x_inv_sqr);
            U = (U_L.pow(&x_sqr) * &U) * U_R.pow(&x_inv_sqr);
            Z = (Z_L.pow(&x_sqr) * &Z) * Z_R.pow(&x_inv_sqr);
            w = w_L;
            v = v_L;
            i += 1;
        }

        //Self::naive_pairing_checks(&proof.A, &proof.B, &w[0], &v[0], &T, &U, &Z)
        Self::batch_pairing_checks(&proof.A, &proof.B, &w[0], &v[0], &T, &U, &Z)
    }
    fn naive_pairing_checks(A: &G1, B: &G2, w: &G1, v: &G2, T: &GT, U: &GT, Z: &GT) -> bool {
        if GT::ate_pairing(A, v) != *T {
            return false
        }
        if GT::ate_pairing(w, B) != *U {
            return false
        }
        if GT::ate_pairing(A, B) != *Z {
            return false
        }
        true
    }

    fn batch_pairing_checks(A: &G1, B: &G2, w: &G1, v: &G2, T: &GT, U: &GT, Z: &GT) -> bool {
        // Check e(A, v) == T and e(w, B) == U and e(A, B) == Z
        // e(A, v) == T
        // e(w, B)^r == U^r
        // e(A, B)^{r^2} == Z^{r^2}
        // e(A, v) * e(w, B)^r * e(A, B)^{r^2} == T * U^r * Z^{r^2}
        // e(A, v) * e(w^r, B) * e(A^{r^2}, B) == T * U^r * Z^{r^2}
        let mut pairs = vec![];
        let r = FieldElement::random();
        let r_sqr = r.square();
        let w_r = w * &r;
        let A_r_sqr = A * &r_sqr;
        pairs.push((A, v));
        pairs.push((&w_r, B));
        pairs.push((&A_r_sqr, B));
        let lhs = GT::ate_multi_pairing(pairs);
        let rhs = &(T * &(U.pow(&r))) * &(Z.pow(&r_sqr));
        lhs == rhs
    }

    pub fn verify_proof(m: usize, transcript: &mut Transcript, T: &GT, U: &GT,
                        Z: &GT, w: &G1Vector, v: &G2Vector, proof: &InnerPairingProductProof) -> bool {
        // TODO:
        // w and v can be computed using a one big multi-scalar multiplication.
        // But not sure of computing T, U and Z from each recursion round in an efficient way since
        // they are members of extension field. This might help https://eprint.iacr.org/2013/458.pdf.
        // Also check pow4 of FP12.rs. It can be generalized.
        unimplemented!()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::time::Instant;

    #[test]
    fn test_inner_pairing_product_argument_recursive_verif() {
        let m = 4;
        let w = G1Vector::random(m);
        let v = G2Vector::random(m);
        let A = G1Vector::random(m);
        let B = G2Vector::random(m);

        let T = GT::inner_product(A.as_slice(), v.as_slice()).unwrap();
        let U = GT::inner_product(w.as_slice(), B.as_slice()).unwrap();
        let Z = GT::inner_product(A.as_slice(), B.as_slice()).unwrap();

        let mut new_trans = Transcript::new(b"innerproduct");
        let start = Instant::now();
        let ipp_proof = InnerPairingProduct::gen_proof(&mut new_trans, &T, &U, &Z, &w, &v, &A, &B);
        println!(
            "Time for create inner product proof for vectors with {} items is {:?}",
            m,
            start.elapsed()
        );

        let mut new_trans = Transcript::new(b"innerproduct");
        let start = Instant::now();
        let res = InnerPairingProduct::verify_proof_recursively(m, &mut new_trans, &T, &U, &Z, &w, &v, &ipp_proof);
        assert!(res);
        println!(
            "Time for verify inner product proof for vectors with {} items is {:?}",
            m,
            start.elapsed()
        );
    }
}