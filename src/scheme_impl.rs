use crate::scheme_traits::SharingScheme;
use crate::hadamard_matrix::HadamardMatrix;
use rand::Rng;
use anyhow::Result;
use ndarray::{arr2, Array2};
use std::mem::size_of;

#[cfg(feature = "zeroize_memory")]
use zeroize::Zeroize;

// что такое Datatype
#[derive(Clone)]
#[derive(Copy)]
#[cfg_attr(feature = "zeroize_memory", derive(Zeroize))]
#[cfg_attr(feature = "zeroize_memory", zeroize(drop))]
pub struct Part {
    number: usize,
    data: u32,
}

impl Part {
    pub fn new(number_: usize, data_: u32) -> Self {
        Part{
            number: number_,
            data: data_,
        }
    }
    pub fn number(&self) -> usize {
        self.number
    }
    pub fn data(&self) -> u32 {
        self.data
    }
}

#[cfg_attr(feature = "zeroize_memory", derive(Zeroize))]
#[cfg_attr(feature = "zeroize_memory", zeroize(drop))]
pub struct HSS {
    mtx: Array2<i32>,
}

impl HSS {
    pub fn new(mtx: &Array2<i32>) -> Self {
        HSS {
            mtx:mtx.clone()
        }
    }

    pub fn mtx_len(&self) -> usize {
        self.mtx.len()
    }

}

impl SharingScheme for HSS {
    type Error = &'static str;
    type SecretType = u32;
    type PartType = Part;

    fn share(&self, secret: Self::SecretType) -> Result<Vec<Part>, Self::Error> {
        let n = self.mtx.shape()[0];
        let mut rng = rand::thread_rng();
        let mut res: Vec<Part> = vec![Part{number: 0, data: 0}; n];
        let secret_size = size_of::<Self::SecretType>() * 8;
        let times = (secret_size + n - 1) / n;
        for i in 0..n {
            res[i].number = i;
            for s_ind in 0..times {
                for j in 0..n {
                    let j_id = j + s_ind * n;
                    if j_id < secret_size {
                        if self.mtx[[i, j]] == 1 {
                            res[i].data |= (1 << j_id) & secret;
                        } else {
                            res[i].data |= (1 << j_id) * rng.gen_range(0..=1);
                        }
                    }
                }
            }
        }
        Ok(res)
    }

    fn reconstruct(&self, parts: Vec<Part>) -> Result<Self::SecretType, Self::Error> {
        let n = self.mtx.shape()[0];
        let mut res: Self::SecretType = 0;
        let secret_size = size_of::<Self::SecretType>() * 8;
        let times = (secret_size + n - 1) / n;
        for i in 0..parts.len() {
            let ind = parts[i].number;
            for s_ind in 0..times {
                for j in 0..n {
                    let j_id = j + s_ind * n;
                    if j_id < secret_size {
                        if self.mtx[[ind, j]] == 1 {
                            res |= (1 << j_id) & parts[i].data;
                        }
                    }
                }
            }
        }
        Ok(res)
    }

    fn validate(&self, parts: Vec<Part>) -> Vec<Part> {
        let n = self.mtx.shape()[0];
        let secret_size = size_of::<Self::SecretType>() * 8;
        let times = (secret_size + n - 1) / n;
        let mut cells: Vec<Vec<Vec<i32>>> = vec![vec![vec![]; 2]; secret_size];
        for i in 0..parts.len() {
            let ind = parts[i].number;
            for s_ind in 0..times {
                for j in 0..n {
                    let j_id = j + s_ind * n;
                    if j_id < secret_size {
                        if self.mtx[[ind, j]] == 1{
                            let bit = (((1 << j_id) & parts[i].data) > 0) as usize;
                            cells[j_id][bit].push(ind as i32);
                        }
                    }
                }
            }
        }

        let mut suspicious: Vec<bool> = vec![false; secret_size];
        for i in 0..secret_size {
            if !cells[i][0].is_empty() && !cells[i][1].is_empty() {
                let more = (cells[i][0].len() > cells[i][1].len()) as usize;
                for ind in &cells[i][more] {
                    suspicious[*ind as usize] = true;
                }
            }
        }

        let mut res: Vec<Part> = Vec::new();
        for i in 0..secret_size {
            if suspicious[i] {
                res.push(parts[i]);
            }
        }
        res
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reconstruction() {
        let h_mtx = HadamardMatrix::from(&arr2(&[[1, 1, 1, 1, 1, 1, 1, 1],
                                                [1, -1, 1, -1, 1, -1, 1, -1],
                                                [1, 1, -1, -1, 1, 1, -1, -1],
                                                [1, -1, -1, 1, 1, -1, -1, 1],
                                                [1, 1, 1, 1, -1, -1, -1, -1],
                                                [1, -1, 1, -1, -1, 1, -1, 1],
                                                [1, 1, -1, -1, -1, -1, 1, 1],
                                                [1, -1, -1, 1, -1, 1, 1, -1]]))
            .unwrap()
            .normalize()
            .get_incidence();
        let hss = HSS::new(&h_mtx);
        for secret in 0..100 {
            let res = hss.share(secret).unwrap();
            let secret_res = hss.reconstruct(res[0..5].to_vec()).unwrap();
            assert_eq!(secret, secret_res);
        }
    }

    #[test]
    fn test_validate() {
        let h_mtx = HadamardMatrix::from(&arr2(&[[1, 1, 1, 1, 1, 1, 1, 1],
                                                [1, -1, 1, -1, 1, -1, 1, -1],
                                                [1, 1, -1, -1, 1, 1, -1, -1],
                                                [1, -1, -1, 1, 1, -1, -1, 1],
                                                [1, 1, 1, 1, -1, -1, -1, -1],
                                                [1, -1, 1, -1, -1, 1, -1, 1],
                                                [1, 1, -1, -1, -1, -1, 1, 1],
                                                [1, -1, -1, 1, -1, 1, 1, -1]]))
            .unwrap()
            .normalize()
            .get_incidence();
        let hss = HSS::new(&h_mtx);
        for secret in 0..100 {
            let res = hss.share(secret).unwrap();
            let valid = hss.validate(res[0..5].to_vec()).is_empty();
            let secret_res = hss.reconstruct(res[0..5].to_vec()).unwrap();
            println!("secret {}, secret_res {}", secret, secret_res);
            assert_eq!(valid, (secret == secret_res));
        }
        for secret in 0..100 {
            let mut res = hss.share(secret).unwrap();
            res[0] = Part::new(res[0].number(), res[0].data() ^ 43);
            let valid = hss.validate(res[0..5].to_vec()).is_empty();
            let secret_res = hss.reconstruct(res[0..5].to_vec()).unwrap();
            assert_eq!(valid, (secret == secret_res));
        }
    }
}
