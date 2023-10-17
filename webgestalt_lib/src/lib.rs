pub mod methods;
pub mod readers;
pub mod stat;
pub enum Error {
    MalformedFile(MalformedError),
    IOError(std::io::Error),
}

pub enum MalformedError {
    NoColumnsFound,
    WrongFormat,
    Unknown,
}

pub enum StatisticsError {
    FoundNANValue,
    InvalidValue,
}

#[cfg(test)]
mod tests {
    use super::*;
}
