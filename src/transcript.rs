use merlin::Transcript;

use amcl_wrapper::field_elem::FieldElement;
use amcl_wrapper::group_elem::GroupElement;
use amcl_wrapper::extension_field_gt::GT;
use amcl_wrapper::constants::FieldElement_SIZE;

pub trait TranscriptProtocol {
    /// Commit a domain separator for a length-`n` inner product proof.
    fn innerproduct_domain_sep(&mut self, n: u64);
    /// Commit a `scalar` with the given `label`.
    fn commit_scalar(&mut self, label: &'static [u8], scalar: &FieldElement);
    /// Commit a GT `point` with the given `label`.
    fn commit_GT(&mut self, label: &'static [u8], point: &GT);
    /// Compute a `label`ed challenge variable.
    fn challenge_scalar(&mut self, label: &'static [u8]) -> FieldElement;
}

impl TranscriptProtocol for Transcript {
    fn innerproduct_domain_sep(&mut self, n: u64) {
        self.append_message(b"dom-sep", b"ipp v1");
        self.append_message(b"n", &n.to_le_bytes());
    }


    fn commit_scalar(&mut self, label: &'static [u8], scalar: &FieldElement) {
        self.append_message(label, &scalar.to_bytes());
    }

    fn commit_GT(&mut self, label: &'static [u8], point: &GT) {
        self.append_message(label, &point.to_bytes());
    }

    fn challenge_scalar(&mut self, label: &'static [u8]) -> FieldElement {
        let mut buf = [0u8; FieldElement_SIZE];
        self.challenge_bytes(label, &mut buf);

        FieldElement::from(&buf)
    }
}