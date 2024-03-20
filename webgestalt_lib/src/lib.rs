use std::{error::Error, fmt};

pub mod methods;
pub mod readers;
pub mod stat;
pub mod writers;

trait CustomError {
    fn msg(&self) -> String;
}

#[derive(Debug)]
pub enum WebGestaltError {
    MalformedFile(MalformedError),
    StatisticsError(StatisticsError),
    IOError(std::io::Error),
}

impl Error for WebGestaltError {}

impl fmt::Display for WebGestaltError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let msg: String = match &self {
            WebGestaltError::MalformedFile(x) => x.msg(),
            WebGestaltError::StatisticsError(x) => x.msg(),
            WebGestaltError::IOError(x) => x.to_string(),
        };
        write!(f, "{}", msg)
    }
}

#[derive(Debug)]
pub struct MalformedError {
    pub path: String,
    pub kind: MalformedErrorType,
}

#[derive(Debug)]
pub enum MalformedErrorType {
    NoColumnsFound { delimeter: String },
    WrongFormat { found: String, expected: String },
    Unknown,
}

impl CustomError for MalformedError {
    fn msg(&self) -> String {
        let error_msg = match &self.kind {
            MalformedErrorType::WrongFormat { found, expected } => format!(
                "Wrong Format Found. Found: {}; Expected: {}",
                found, expected
            ),
            MalformedErrorType::Unknown => String::from("Unknown error type."),
            MalformedErrorType::NoColumnsFound { delimeter } => format!(
                "No column found with delimeter {}",
                if delimeter == "\t" { "\\t" } else { delimeter }
            ),
        };
        format!("Error in {}: {}.", self.path, error_msg)
    }
}

#[derive(Debug)]
pub enum StatisticsError {
    FoundNANValue,
    InvalidValue { value: f64 },
}

impl CustomError for StatisticsError {
    fn msg(&self) -> String {
        let error_msg = match &self {
            StatisticsError::FoundNANValue => String::from("Found a NAN value"),
            StatisticsError::InvalidValue { value } => format!("Found invalid value: {}", value),
        };
        format!("Statstical Error: {}.", error_msg)
    }
}
