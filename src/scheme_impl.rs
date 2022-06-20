//! Реализация библиотечных структур.

use crate::scheme_traits::SharingScheme;
use crate::hadamard_matrix::HadamardMatrix;
use rand::Rng;
use anyhow::Result;
use ndarray::{arr2, Array2};
use std::mem::size_of;

#[cfg(feature = "zeroize_memory")]
use zeroize::Zeroize;

#[derive(Clone)]
#[derive(Copy)]
#[cfg_attr(feature = "zeroize_memory", derive(Zeroize))]
#[cfg_attr(feature = "zeroize_memory", zeroize(drop))]
/// Структура отдельной доли, получаемой при разделении секрета.
pub struct Part {
    /// Номер, соответвуюший строке матрицы Адамара, по которой была получена доля.
    number: usize,
    /// Значение доли.
    data: u32,
}

/// Реализация методов структуры доли.
impl Part {
    /// Создание экземпляра структуры [scheme_impl::Part] по данному номеру и значению.
    pub fn from(number_: usize, data_: u32) -> Self {
        Part{
            number: number_,
            data: data_,
        }
    }

    /// Возвращение значения поля number.
    pub fn number(&self) -> usize {
        self.number
    }

    /// Возвращение значения поля data.
    pub fn data(&self) -> u32 {
        self.data
    }
}

#[cfg_attr(feature = "zeroize_memory", derive(Zeroize))]
#[cfg_attr(feature = "zeroize_memory", zeroize(drop))]
/// Структура схемы разделения секрета, содержащая поле с матрицей инцидентности.
pub struct HSS {
    /// Матрица инцидентности, построенная по матрице Адамара.
    mtx: Array2<i32>,
}

/// Реализация базовых методов структуры схемы разделения секрета.
impl HSS {
    /// Создание экземпляра структуры по данной матрице инцидентности.
    pub fn from(mtx: &Array2<i32>) -> Self {
        HSS {
            mtx:mtx.clone()
        }
    }

    /// Возвращение размерности хранимой матрицы -- максимального числа долей, на которые будет разбит секрет.
    pub fn mtx_len(&self) -> usize {
        self.mtx.len()
    }
}

/// Реализация методов трейта [share_traits::SharingScheme] в структуре [share_impl::HSS].
impl SharingScheme for HSS {
    type Error = &'static str;
    type SecretType = u32;
    type PartType = Part;

    /// Метод, реализующий разбиение секрета типа u32 на n долей типа [share_impl::Part].
    ///
    /// В цикле по i обрабатывается i-я строка матрицы инцидентности.
    ///
    /// Цикл по s_ind обеспечивает пробегание по всем 32-м битам секрета, в случае размерности матрицы инцидентности меньшей 32.
    ///
    /// Рассмотрим, что происходит с j_id-м битом секрета (j_id = j + s_ind * n):
    /// - mtx[[i, j]] == 1, j_id-й бит приравнивается j_id-му биту секрета
    /// - mtx[[i, j]] == 0, j_id-й бит приравнивается рандомному значению {0, 1}
    ///
    /// # Пример.
    /// ```
    /// let res = hss.share(secret).unwrap();
    /// ```
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

    /// Восстановление секрета по данному набору долей. Не происходит никаких проверок. Как следствие, в случае ошибки в какой-то доли, восстановленный секрет может отличаться от исходного.
    ///
    /// Проходимся по строке матрицы Адамара, если в j-м элементе стоит 1, то в итоговом значении секрета соответствующему j_id-му биту проставляем j_id-й бит из доли.
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

    /// Проверка на корректность пришедшего набора долей.
    ///
    /// Формируем трёхмерный вектор cells[bit_number][bit_value][part_number] хранящий
    /// информацию для каждого j_id-го бита, номера каких частей дают значение 1, а каких 0.
    ///
    /// Далее по вектору cells определяем подозрительные части -- в вектор флагов suspicious
    /// в случае присутствия одновременно номеров долей в cells[i][0] и cells[i][1]
    /// проставляем 1 в индексах элементов, соответвуюших долям из наименьшего множества
    /// (из cells[i][0] или cells[i][1]).
    ///
    /// По проставленным флагам формируем вектор, хранящий номера подозрительных долей.
    fn validate(&self, parts: Vec<Part>) -> Vec<usize> {
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

        let mut res: Vec<usize> = Vec::new();
        for i in 0..secret_size {
            if suspicious[i] {
                res.push(i);
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
        let hss = HSS::from(&h_mtx);
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
        let hss = HSS::from(&h_mtx);
        for secret in 0..100 {
            let res = hss.share(secret).unwrap();
            let valid = hss.validate(res[0..5].to_vec()).is_empty();
            let secret_res = hss.reconstruct(res[0..5].to_vec()).unwrap();
            println!("secret {}, secret_res {}", secret, secret_res);
            assert_eq!(valid, (secret == secret_res));
        }
        for secret in 0..100 {
            let mut res = hss.share(secret).unwrap();
            res[0] = Part::from(res[0].number(), res[0].data() ^ 43);
            let valid = hss.validate(res[0..5].to_vec()).is_empty();
            let secret_res = hss.reconstruct(res[0..5].to_vec()).unwrap();
            assert_eq!(valid, (secret == secret_res));
        }
    }
}
