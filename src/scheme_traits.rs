//! Трейты, необходимые для реализации в библиотеке.

use std::vec::Vec;
use anyhow::Result;

pub trait SharingScheme {
    type Error;
    type SecretType;
    type PartType;
    /// Разделение секрета на доли.
    fn share(&self, secret: Self::SecretType) -> Result<Vec<Self::PartType>, Self::Error>;
    /// Восстановление секрета по вектору долей.
    fn reconstruct(&self, shares: Vec<Self::PartType>) -> Result<Self::SecretType, Self::Error>;
    /// Валидация множества долей: возвращение номеров долей, предположительно
    /// используемых злоумышленниками.
    fn validate(&self, shares: Vec<Self::PartType>) -> Vec<usize>;
}
