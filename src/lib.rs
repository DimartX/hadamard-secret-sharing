//! # Hadamard secret sharing scheme
//!
//! `hadamard_sss` это библиотека, реализующая схему разделения секрета
//! на основе матриц Адамара.
//!
//! Если кратко, реализованы методы трейта [scheme_traits::SharingScheme]:
//! ```
//! fn share(&self, secret: Self::SecretType) -> Result<Vec<Self::PartType>, Self::Error>;
//! fn reconstruct(&self, shares: Vec<Self::PartType>) -> Result<Self::SecretType, Self::Error>;
//! fn validate(&self, shares: Vec<Self::PartType>) -> Vec<Self::PartType>;
//! ```
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

/// Основная структура
pub struct HadamardSSS {
    /// внутренняя структура реализующая схему разделения секрета
    hss: HSS,
    /// пороговое значение для матрицы Адамара, переданной в структуру
    threshold: usize,
}

/// Реализация структуры
impl HadamardSSS {
    /// Создание экземпляра структуры по данной матрице Адамара
    pub fn from(mtx: &Array2<i32>) -> Result<Self, &'static str> {
        let mut had = HadamardMatrix::from(&mtx).expect("Error! ");
        let incidence_mtx = had.normalize().get_incidence();
        Ok(HadamardSSS {
            hss: HSS::from(&incidence_mtx),
            threshold: HadamardSSS::get_threshold(&incidence_mtx),
        })
    }

    /// Возвращение порогового числа участников, необходимого для восстановления секрета
    pub fn get_threshold(mtx: &Array2<i32>) -> usize {
        // Соображения: mtx.shape()[0] == 4n - 1, threshold = 2n + 1 = (4n - 1 + 3) / 2
        (mtx.shape()[0] + 3) / 2
    }

    /// Проверка, образуют ли представленные доли валидный набор для восстановления секрета
    pub fn is_valid(&self, parts: Vec<Part>) -> bool {
        self.hss.validate(parts).is_empty()
    }
}

/// Реализация трейта SharingScheme в структуре HadamardSSS
impl SharingScheme for HadamardSSS {
    type Error = &'static str;
    type SecretType = u32;
    type PartType = Part;

    /// Обёртка для share_impl::HSS::share
    fn share(&self, secret: Self::SecretType) -> Result<Vec<Part>, Self::Error> {
        self.hss.share(secret)
    }

    /// Обёртка для share_impl::HSS::reconstruct с учётом количества пришёдших долей
    fn reconstruct(&self, parts: Vec<Self::PartType>) -> Result<Self::SecretType, &'static str> {
        if parts.len() < self.threshold {
            println!("{} is less than threshold {} parties", parts.len(), self.threshold);
            Err("less than threshold parties")
        } else {
            self.hss.reconstruct(parts)
        }
    }

    /// Обёртка для share_impl::HSS::validate
    fn validate(&self, parts: Vec<Part>) -> Vec<usize> {
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
        let hsss = HadamardSSS::from(&h_mtx).unwrap();
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
        let hsss = HadamardSSS::from(&h_mtx).unwrap();
        for secret in 0..100 {
            let res = hsss.share(secret).unwrap();
            let valid = hsss.validate(res[0..5].to_vec()).is_empty();
            let secret_res = hsss.reconstruct(res[0..5].to_vec()).unwrap();
            assert_eq!(valid, (secret == secret_res));
        }
        for secret in 0..100 {
            let mut res = hsss.share(secret).unwrap();
            res[0] = Part::from(res[0].number(), res[0].data() ^ 43);
            let valid = hsss.validate(res[0..5].to_vec()).is_empty();
            let secret_res = hsss.reconstruct(res[0..5].to_vec()).unwrap();
            assert_eq!(valid, (secret == secret_res));
        }
    }
}
