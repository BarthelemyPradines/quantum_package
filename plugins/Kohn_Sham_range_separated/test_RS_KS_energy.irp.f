program test_RS_KS
 implicit none
 read_wf = .False.
 touch read_wf

 call test_energy 
end

subroutine test_energy
 implicit none
 print*, 'exchange_functional    =  ',exchange_functional
 print*, 'correlation_functional =  ',correlation_functional
 print*, 'SCF_energy          = ',SCF_energy
 print*, 'RS_KS_energy        = ',RS_KS_energy
 print*, 'one_electron_energy = ',one_electron_energy
 print*, 'two_electron_energy = ',two_electron_energy
 print*, 'e_exchange_dft      = ',e_exchange_dft
 print*, 'e_correlation_dft   = ',e_correlation_dft
 print*, 'Fock_matrix_energy  = ',Fock_matrix_energy
end


