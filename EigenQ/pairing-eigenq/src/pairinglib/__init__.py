"""pairinglib: resolution-refinement + rodeo pipeline for the pairing Hamiltonian."""
from .hamiltonian import build_sector, H_pairing_sparse, fci_ground, E_HF, build_H_full
from .ccd import run_ccd
from .uccsd import uccsd_pool, setup_uccsd, uccsd_vqe
from .gates import apply_1q, apply_cnot, Ry, Rz, pauli_exp
from .gate_uccsd import gate_uccsd_setup, gate_uccsd_vqe, gate_uccsd_state
from .refine import prolong_matrix, refine_state, refine_NK, refine_from_state
from .rodeo import rodeo_post0, rodeo_cycle_ancilla, rodeo_sweep, rodeo_track, rodeo_allzero_prob
from .trotter import (pairing_terms, ph_terms, model_terms, trotter_U, trotter_U_dt,
                      trotter_error, refine_state_trotter, rodeo_track_U, cnot_counts,
                      ph_pauli_costs)
from .phmodel import (H_ph_sparse, build_H_ph_full, E_ref, seniority_of,
                      seniority_weight, seniority_zero_mask, ph_coupling_table,
                      ph_term_sector, ph_term_full)
from .pauli import pauli_decompose, pauli_count

__version__ = "0.1.0"
