//! Модуль, в котором реализована структура для работы с матрицами Адамара.
use ndarray::{arr2, Array2};
#[cfg(feature = "zeroize_memory")]
use zeroize::Zeroize;

#[cfg_attr(feature = "zeroize_memory", derive(Zeroize))]
#[cfg_attr(feature = "zeroize_memory", zeroize(drop))]
/// Структура, хранящая матрицу Адамара.
pub struct HadamardMatrix {
    /// Двумерная матрица, являющаяся матрицей Адамара.
    mtx: Array2<i32>,
}

/// Реализация методов структуры.
impl HadamardMatrix {
    /// Создание экземпляра структуры по данной двумерной матрице.
    /// В случае, если пришедшая матрица не является матрицей Адамара,
    /// возвращается соответствующая ошибка.
    /// # Пример
    /// ```
    /// let h_mtx = HadamardMatrix::from(&arr2(&[[1, 1],
    ///                                         [1, -1]]).expect("Can't create Hadamard mtx."));
    /// ```
    pub fn from(mtx: &Array2<i32>) -> Result<HadamardMatrix, &'static str> {
        if HadamardMatrix::is_hadamard(&mtx) {
            Ok(HadamardMatrix { mtx: mtx.clone() })
        } else {
            Err("something wrong with that matrix")
        }
    }

    /// Проверка матрицу на Адамаровость:
    /// - является ли она квадратной
    /// - является ли она непустой
    /// - состоит ли только из -1 и 1
    /// - проверка на определение H * H.T = nI
    /// # Пример.
    /// ```
    /// let res = HadamardMatrix::is_hadamard(&arr2(&[[1, 2],
    ///                                               [3, 4]]);
    /// // получим res == false
    fn is_hadamard(mtx: &Array2<i32>) -> bool {
        let n = mtx.shape()[0];
        let m = mtx.shape()[1];
        if !mtx.is_square() || n < 1 {
            return false;
        }
        for i in 0..n {
            for j in 0..n {
                if (mtx[[i, j]] != -1) && (mtx[[i, j]] != 1) {
                    return false;
                }
            }
        }
        let res = n as i32 * Array2::<i32>::eye(n);
        let mult_res = mtx.dot(&mtx.t());

        if res != mult_res {
            false
        } else {
            true
        }
    }

    /// Нормализация матрицы Адамара, чтобы первый столбец и первая строка состояли из одних 1.
    /// ```
    /// [[-1, -1],
    ///  [-1, 1]]
    /// ```
    /// Становится
    /// ```
    /// [[1, 1],
    ///  [1, -1]]
    /// ```
    pub fn normalize(&mut self) -> &mut Self {
        let n = self.mtx.shape()[0];
        for i in 0..n {
            if self.mtx[[i, 0]] == -1 {
                for j in 0..n {
                    self.mtx[[i, j]] *= -1;
                }
            }
            if self.mtx[[0, i]] == -1 {
                for j in 0..n {
                    self.mtx[[j, i]] *= -1;
                }
            }
        }
        self
    }

    /// Получения матрицы инцидентности по данной матрице Адамара,
    /// которая будет соответствовать блок-дизайну 2-(4n-1, 2n-1, n-1)
    pub fn get_incidence(&self) -> Array2<i32> {
        let n = self.mtx.shape()[0];
        (&self.mtx.view().slice(s![1.., 1..]) + &Array2::<i32>::ones((n - 1, n - 1))) / 2
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic]
    fn empty_mtx() {
        HadamardMatrix::from(&arr2(&[[]])).expect("this is not hadamard matrix");
    }

    #[test]
    fn false_mtx() {
        assert_eq!(HadamardMatrix::is_hadamard(&arr2(&[[1, 2],
                                                   [3, 4]]))
                   , false);
        assert_eq!(HadamardMatrix::is_hadamard(&arr2(&[[1, 0],
                                                    [0, 1]]))
                   , false);
        assert_eq!(HadamardMatrix::is_hadamard(&arr2(&[[-1, -1],
                                                    [1, 1]]))
                   , false);
    }

    #[test]
    fn true_mtx() {
        assert_eq!( HadamardMatrix::is_hadamard(&arr2(&[[1, 1],
                                                   [1, -1]]))
                   , true);
        assert_eq!(HadamardMatrix::is_hadamard(&arr2(&[[1, 1, 1, 1],
                                                    [1, -1, 1, -1],
                                                    [1, 1, -1, -1],
                                                    [1, -1, -1, 1]]))
                   , true);
    }

    #[test]
    fn test_normalize() {
        assert_eq!(HadamardMatrix::from(&arr2(&[[-1, -1],
                                               [-1, 1]])).unwrap().normalize().mtx
                   ,
                   HadamardMatrix::from(&arr2(&[[1, 1],
                                               [1, -1]])).unwrap().mtx
        );
        assert_eq!(HadamardMatrix::from(&arr2(&[[-1, 1],
                                               [1, 1]])).unwrap().normalize().mtx
                   ,
                   HadamardMatrix::from(&arr2(&[[1, 1],
                                               [1, -1]])).unwrap().mtx
        );
    }


    #[test]
    fn test_incidence() {
        assert_eq!(HadamardMatrix::from(&arr2(&[[1, 1],
                                               [1, -1]])).unwrap().get_incidence()
                   ,
                   arr2(&[[0]])
        );
        assert_eq!(HadamardMatrix::from(&arr2(&[[1, 1, 1, 1],
                                               [1, -1, 1, -1],
                                               [1, 1, -1, -1],
                                               [1, -1, -1, 1]])).unwrap().get_incidence()
                   ,
                   arr2(&[[0, 1, 0],
                          [1, 0, 0],
                          [0, 0, 1]])
        );
    }
}
