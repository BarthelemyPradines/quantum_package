program pouet
 implicit none
 read_wf = .true.
 touch read_wf
 call routine2

end
subroutine routine2
 implicit none

 print*,  '****************************************'
 write(*, '(A22,X,F16.10)') 'TOTAL ENERGY        = ',total_electronic_energy
 print*, ''
 print*, 'Component of the energy ....'
 print*, ''
 write(*, '(A22,X,F16.10)') 'psi_energy_erf      = ',psi_energy_erf      
 write(*, '(A22,X,F16.10)') 'psi_energy_core     = ',psi_energy_core    
 write(*, '(A22,X,F16.10)') 'psi_energy_hartree  = ',psi_energy_hartree
 write(*, '(A22,X,F16.10)') 'energy_x            = ',energy_x         
 write(*, '(A22,X,F16.10)') 'energy_c            = ',energy_c          
 print*, ''
 print*,  '****************************************'
 print*, ''
 write(*, '(A22,X,F16.10)') 'Approx eigenvalue   = ',Fock_matrix_expectation_value + psi_energy_erf
 write(*, '(A22,X,F16.10)') 'Trace_v_xc          = ',Trace_v_xc


end
