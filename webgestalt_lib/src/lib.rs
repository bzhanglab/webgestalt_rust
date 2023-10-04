pub mod methods;
pub mod readers;

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
}

#[cfg(test)]
mod tests {
    use super::*;
}
