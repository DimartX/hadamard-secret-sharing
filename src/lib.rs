#[macro_use]
extern crate ndarray;
mod scheme_impl;
mod scheme_traits;
mod hadamard_matrix;
use hadamard_matrix::HadamardMatrix;
use scheme_impl::{HSS, Part};
use scheme_traits::SharingScheme;
use anyhow::Result;
use ndarray::{arr2, Array2};

struct HadamardSSS {
    hss: HSS,
    threshold: usize,
}

impl HadamardSSS {
    pub fn new(mtx: &Array2<i32>) -> Result<Self, &'static str> {
        let mut had = HadamardMatrix::from(&mtx).expect("Error! ");
        let incidence_mtx = had.normalize().get_incidence();
        Ok(HadamardSSS {
            hss: HSS::new(&incidence_mtx),
            threshold: HadamardSSS::get_threshold(&incidence_mtx),
        })
    }

    fn get_threshold(mtx: &Array2<i32>) -> usize {
        // this: mtx.shape()[0] == 4n - 1, threshold = 2n + 1 = (4n - 1 + 3) / 2
        (mtx.shape()[0] + 3) / 2
    }

    fn is_valid(&self, parts: Vec<Part>) -> bool {
        self.hss.validate(parts).is_empty()
    }

}

impl SharingScheme for HadamardSSS {
    type Error = &'static str;
    type SecretType = u32;
    type PartType = Part;

    fn share(&self, secret: Self::SecretType) -> Result<Vec<Part>, Self::Error> {
        self.hss.share(secret)
    }

    fn reconstruct(&self, parts: Vec<Self::PartType>) -> Result<Self::SecretType, &'static str> {
        if parts.len() < self.threshold {
            println!("{} is less than threshold {} parties", parts.len(), self.threshold);
            Err("less than threshold parties")
        } else {
            self.hss.reconstruct(parts)
        }
    }

    fn validate(&self, parts: Vec<Part>) -> Vec<Part> {
        self.hss.validate(parts)
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_secret_reconstruction() {
        let h_mtx = arr2(&[[1, 1, 1, 1, 1, 1, 1, 1],
                           [1, -1, 1, -1, 1, -1, 1, -1],
                           [1, 1, -1, -1, 1, 1, -1, -1],
                           [1, -1, -1, 1, 1, -1, -1, 1],
                           [1, 1, 1, 1, -1, -1, -1, -1],
                           [1, -1, 1, -1, -1, 1, -1, 1],
                           [1, 1, -1, -1, -1, -1, 1, 1],
                           [1, -1, -1, 1, -1, 1, 1, -1]]);
        let hsss = HadamardSSS::new(&h_mtx).unwrap();
        for secret in 0..100 {
            let res = hsss.share(secret).unwrap();
            let secret_res = hsss.reconstruct(res[0..5].to_vec()).unwrap();
            assert_eq!(secret, secret_res);
        }
    }

    #[test]
    fn test_hsss_validate() {
        let h_mtx = arr2(&[[1, 1, 1, 1, 1, 1, 1, 1],
                           [1, -1, 1, -1, 1, -1, 1, -1],
                           [1, 1, -1, -1, 1, 1, -1, -1],
                           [1, -1, -1, 1, 1, -1, -1, 1],
                           [1, 1, 1, 1, -1, -1, -1, -1],
                           [1, -1, 1, -1, -1, 1, -1, 1],
                           [1, 1, -1, -1, -1, -1, 1, 1],
                           [1, -1, -1, 1, -1, 1, 1, -1]]);
        let hsss = HadamardSSS::new(&h_mtx).unwrap();
        for secret in 0..100 {
            let res = hsss.share(secret).unwrap();
            let valid = hsss.validate(res[0..5].to_vec()).is_empty();
            let secret_res = hsss.reconstruct(res[0..5].to_vec()).unwrap();
            assert_eq!(valid, (secret == secret_res));
        }
        for secret in 0..100 {
            let mut res = hsss.share(secret).unwrap();
            res[0] = Part::new(res[0].number(), res[0].data() ^ 43);
            let valid = hsss.validate(res[0..5].to_vec()).is_empty();
            let secret_res = hsss.reconstruct(res[0..5].to_vec()).unwrap();
            assert_eq!(valid, (secret == secret_res));
        }
    }
}
