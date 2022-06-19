use std::vec::Vec;
use anyhow::Result;

pub trait SharingScheme {
    type Error;
    type SecretType;
    type PartType;
    fn share(&self, secret: Self::SecretType) -> Result<Vec<Self::PartType>, Self::Error>;
    fn reconstruct(&self, shares: Vec<Self::PartType>) -> Result<Self::SecretType, Self::Error>;
    fn validate(&self, shares: Vec<Self::PartType>) -> Vec<Self::PartType>;
}
