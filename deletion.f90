!********************************************************************************
!   Cassandra - An open source atomistic Monte Carlo software package
!   developed at the University of Notre Dame.
!   http://cassandra.nd.edu
!   Prof. Edward Maginn <ed@nd.edu>
!   Copyright (2013) University of Notre Dame du Lac
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*******************************************************************************
SUBROUTINE Deletion(this_box,mcstep,randno)
  
  !********************************************************************************
  ! This subroutine deletes a molecule a box. Thexs box is first chosen at random
  ! and then an adsorbate component is selected at random. Next a randomly chosen
  ! molecule of this component is deleted from the box
  !
  ! Called by
  !
  !  gcmc_driver
  ! 
  ! Revision History
  !
  !   12/10/13 : Beta Release
  !
  !
  !********************************************************************************
  
  USE Type_Definitions
  USE Run_Variables
  USE Random_Generators
  USE Simulation_Properties
  USE Energy_Routines
  USE Fragment_Growth

  IMPLICIT NONE


  INTEGER, INTENT(INOUT) :: this_box

  INTEGER :: is, k, position, i_type, i, nsorbate, nmols_sorbate, which_anchor
  INTEGER :: is_1, is_counter, mcstep, ktothen, ifrag
  INTEGER, ALLOCATABLE :: sorbate_id(:), frag_order(:)

  REAL(DP) :: delta_e, delta_e_pacc, delta_v, E_bond, E_angle, E_dihedral, E_intra_vdw, E_intra_qq, E_inter_vdw
  REAL(DP) :: E_inter_qq, E_improper
  REAL(DP) :: E_reciprocal_move, E_self_move, e_lrc, alpha_ratio, factor
  REAL(DP) :: gnew, gold, dg ! weighting factors for the new mol number and old mol number
  REAL(DP), ALLOCATABLE :: sorbate_x(:)
  REAL(DP) :: pick_species, this_lambda, randno
  REAL(DP) :: E_intra_vdw_igas, E_intra_qq_igas

  LOGICAL :: accept, accept_or_reject


  REAL(DP) :: check_e, energy_old
  LOGICAL  :: superbad
!-----------------------------------------------------------------------------------  
!FSL
  INTEGER :: alive(nspecies), sm, am, im(nspecies)
  REAL(DP) :: P_reverse(nspecies), nrg_ring_frag_tot(nspecies)
  LOGICAL :: cbmc_overlap(nspecies), intra_overlap(nspecies)
  LOGICAL :: inter_overlap(nspecies), p_overlap
  REAL(DP) :: tb, ta, td, ti, tinterqq, tintervdw, tintraqq, tintravdw, tself
  REAL(DP) :: p_vdw, p_qq, tpacc, ntemp, suben, adden

!FSL
!-----------------------------------------------------------------------------------
! Solid framework

  !   End of section

!RETURN
  
  pacc = 0.0_DP
  paccbiased = 0.0_DP
  alpha_ratio = 1.0_DP
  nrg_ring_frag_tot(:) = 0.0_DP
  im(:) = 0.0_DP
  am = 0.0_DP
  sm = 0.0_DP
  is = 0
  

  ! Randomly choose a box to delete particle from

  IF (int_sim_type .NE. sim_gemc) THEN
     this_box = INT( nbr_boxes * rranf()) + 1
     ! else the box is specified as input in a GEMC simulation
  END IF

  tot_trials(this_box) = tot_trials(this_box) + 1

!  is_1 = INT(rranf() * nspec_insert) + 1
!  is_counter = 0

!  DO is = 1, nspecies
!     IF(species_list(is)%int_species_type == int_sorbate) is_counter = is_counter + 1
!     IF(is_counter == is_1) EXIT
!  END DO

!!$  IF(int_sim_type == sim_gcmc) THEN
!!$     IF( is == 1 ) t_collect = .TRUE.
!!$  END IF

  do is = 1, nspecies
  IF(nmols(is,this_box) == 0) RETURN
  enddo

  ntrials(is,this_box)%deletion = ntrials(is,this_box)%deletion + 1  
  ! Determine the index of im

  do is = 1, nspecies
  im(is) = INT(rranf() * nmols(is,this_box)) + 1

  CALL Get_Index_Molecule(this_box,is,im(is),alive(is))
  if(is == 1) sm = alive(is)
  if(is == 2) am = alive(is)
  enddo
  
  ! Compute the energy of the molecule
  
  ! Intra molecule energy
  
  delta_e = 0.0_DP
  P_reverse(:) = 1.0_DP

  this_lambda = 1.0_DP

  cbmc_overlap(:) = .FALSE.

  do is = 1, nspecies
  IF(species_list(is)%fragment .AND. (species_list(is)%int_insert .NE. int_igas)) THEN
     del_Flag = .TRUE.
     get_fragorder = .TRUE.
     ALLOCATE(frag_order(nfragments(is)))
     CALL Build_Molecule(alive(is),is,this_box,frag_order,this_lambda, which_anchor, P_reverse(is), nrg_ring_frag_tot(is), cbmc_overlap(is))
     DEALLOCATE(frag_order)
     
     IF (cbmc_overlap(is)) THEN
        CALL Revert_Old_Cartesian_Coordinates(alive(is),is)
        atom_list(1:natoms(is),alive(is),is)%exist = .TRUE.
        molecule_list(alive(is),is)%cfc_lambda = this_lambda
        
        WRITE(*,*)
        WRITE(*,*) 'Warning....energy overlap detected in old configuration in deletion.f90'
        WRITE(*,*) 'molecule, species', alive(is), is
        WRITE(*,*)
     END IF

     ktothen = 1

     IF (nfragments(is) /=0 ) THEN

        ktothen = ktothen * kappa_ins

        IF (kappa_rot /= 0) THEN
           ktothen = ktothen * kappa_rot
        END IF
        
        IF (kappa_dih /=0 ) THEN

           DO ifrag = 2, nfragments(is)
              ktothen = ktothen * kappa_dih
           END DO
        END IF
     
     END IF
     
     P_reverse(is) = P_reverse(is) * REAL(ktothen , DP)

  END IF
  enddo

  E_bond = 0.0_DP
  E_angle = 0.0_DP
  E_dihedral = 0.0_DP
  E_improper = 0.0_DP
  do is = 1, nspecies
  tb = 0.0_DP
  ta = 0.0_DP
  td = 0.0_DP
  ti = 0.0_DP
  CALL Get_COM(alive(is),is)
  CALL Compute_Max_COM_Distance(alive(is),is)

  CALL Compute_Molecule_Bond_Energy(alive(is),is,tb)
  CALL Compute_Molecule_Angle_Energy(alive(is),is,ta)
  CALL Compute_Molecule_Dihedral_Energy(alive(is),is,td)
  CALL Compute_Molecule_Improper_Energy(alive(is),is,ti)

  E_bond = E_bond + tb
  E_angle = E_angle + ta
  E_dihedral = E_dihedral + td
  E_improper = E_improper + ti
  enddo

  delta_e = delta_e + E_bond + E_angle + E_dihedral + E_improper

  E_intra_qq = 0.0_DP
  E_intra_vdw = 0.0_DP 
  E_inter_qq = 0.0_DP
  E_inter_vdw = 0.0_DP
  do is = 1, nspecies
  tinterqq = 0.0_DP
  tintervdw = 0.0_DP
  tintraqq = 0.0_DP
  tintravdw = 0.0_DP
  ! Nonbonded energy  
  CALL Compute_Molecule_Nonbond_Intra_Energy(alive(is),is,tintravdw,tintraqq,intra_overlap(is))

!  IF (l_pair_nrg) THEN
!     CALL Store_Molecule_Pair_Interaction_Arrays(alive(is),is,this_box,tintervdw,tinterqq)
!     DEALLOCATE(pair_vdw_temp,pair_qq_temp)
!  ELSE
     CALL Compute_Molecule_Nonbond_Inter_Energy(alive(is),is,tintervdw,tinterqq,inter_overlap(is))
!  END IF
  E_intra_vdw = E_intra_vdw + tintravdw
  E_intra_qq = E_intra_qq + tintraqq
  E_inter_vdw = E_inter_vdw + tintervdw
  E_inter_qq = E_inter_qq + tinterqq
  enddo

!FSL Subtract pair energy due to double counting

       CALL Compute_Molecule_Pair_Interaction(sm,1,am,2,this_box, &
                  p_vdw,p_qq,p_overlap)

!FSL Fim sub pair en double count

  E_inter_vdw = E_inter_vdw - p_vdw
  E_inter_qq = E_inter_qq - p_qq

  delta_e = delta_e + E_intra_vdw + E_intra_qq
  delta_e = delta_e + E_inter_vdw + E_inter_qq

  E_reciprocal_move = 0.0_DP
  E_self_move = 0.0_DP
  do is = 1, nspecies
  tself = 0.0_DP
  IF ( (int_charge_sum_style(this_box) == charge_ewald) .AND. (has_charge(is)) ) THEN
     CALL Ins_Pairs_Ewald_Reciprocal_Energy_Difference(alive(is),alive(is),is,this_box,int_deletion,E_reciprocal_move)
     CALL Compute_Ewald_Self_Energy_Difference(alive(is),is,this_box,int_deletion,tself)
     E_self_move = E_self_move + tself
  END IF
  enddo

     delta_e = delta_e - E_self_move
     delta_e = delta_e - (E_reciprocal_move - energy(this_box)%ewald_reciprocal)

  IF (int_vdw_sum_style(this_box) == vdw_cut_tail) THEN

     ! subtract off beads for this species

     nbeads_out(:) = nint_beads(:,this_box)

     do is = 1, nspecies
     DO i = 1, natoms(is)
        i_type = nonbond_list(i,is)%atom_type_number
        nint_beads(i_type,this_box) = nint_beads(i_type,this_box) - 1
     END DO
     enddo

     CALL Compute_LR_correction(this_box,e_lrc)
     delta_e = delta_e - ( e_lrc - energy(this_box)%lrc )

  END IF  
  

  ! calculate the factor in the acceptance rule and the weighting function in necessary
  ! min(1,exp(-factor)). Note that the energy difference is actually (e_total_new - e_total_old)
  ! delta_e computed above is (e_total_old - e_total_new) so in the factor below
  ! - delta_e is used to represent actual change in energy. See below when energies are
  ! update upon suceessful deletion


    dg = 0.0_DP

    ntemp = 0.0_DP
    suben = 0.0_DP
    adden = 0.0_DP
!FSL energy    
  do is = 1, nspecies
  IF(species_list(is)%int_insert == int_igas) THEN
     igas_flag = .TRUE.
     CALL Compute_Molecule_Nonbond_Intra_Energy(alive(is),is,E_intra_vdw_igas,E_intra_qq_igas,intra_overlap(is))
     igas_flag = .FALSE.
     tpacc = beta(this_box)*(E_intra_vdw_igas + E_intra_qq_igas)
  ELSEIF (species_list(is)%fragment) THEN
     tpacc = beta(this_box) * (nrg_ring_frag_tot(is))
  ELSE
     tpacc = 0.0_DP
  END IF
  if(is == 1) ntemp = nmols(is,this_box) + 1
  if(is == 1) adden = adden + DLOG(P_reverse(1)+P_reverse(2))
  suben = suben + tpacc + beta(this_box)*E_angle*0.5_DP
  enddo
!FSL

  pacc = beta(this_box)*(-delta_e) - suben - adden - DLOG(ntemp)

  pacc = pacc - DLOG(alpha_ratio) + DLOG(box_list(this_box)%volume)

  is = 1
  IF(lchempot) THEN
     pacc = pacc + beta(this_box) * species_list(1)%chem_potential 
     pacc = pacc - 3.0_DP * DLOG(dbpair)
     pacc = pacc + DLOG(species_list(is)%zig_by_omega)
  ELSE
     pacc = pacc + DLOG(species_list(is)%fugacity * beta(this_box)) - &
                   DLOG(species_list(is)%zig_by_omega)
  END IF
!  write(59,*), E_inter_vdw, E_inter_qq, E_bond, E_angle, E_dihedral, &
!     E_improper, E_intra_vdw, E_intra_qq, E_self_move, E_reciprocal_move, &
!     -energy(this_box)%ewald_reciprocal, e_lrc, -energy(this_box)%lrc, &
!     -p_vdw, -p_qq
  ! correct for ring biasing if any

  factor = pacc + dg    
!  factor = pacc 
!  write(*,*) 'deletion factor', factor

  accept = accept_or_reject(factor)

  factor = -factor
  IF( factor < 0.0_DP ) THEN
    factor = DEXP(factor)
  ELSE
    factor = 1.0_DP
  END IF

  paccbiased = factor

  IF(cbmc_overlap(1) .OR. cbmc_overlap(2)) accept = .FALSE.

  IF (accept) THEN
     ! Update energies

     energy(this_box)%total = energy(this_box)%total - delta_e
     energy(this_box)%intra = energy(this_box)%intra - E_bond - E_angle - E_dihedral
     energy(this_box)%bond = energy(this_box)%bond - E_bond
     energy(this_box)%angle = energy(this_box)%angle - E_angle
     energy(this_box)%dihedral = energy(this_box)%dihedral - E_dihedral
     energy(this_box)%improper = energy(this_box)%improper - E_improper
     energy(this_box)%intra_vdw = energy(this_box)%intra_vdw - E_intra_vdw
     energy(this_box)%intra_q = energy(this_box)%intra_q - E_intra_qq
     energy(this_box)%inter_vdw = energy(this_box)%inter_vdw - E_inter_vdw
     energy(this_box)%inter_q   = energy(this_box)%inter_q - E_inter_qq

     is = 1
     IF ( int_charge_sum_style(this_box) == charge_ewald .AND. has_charge(is)) THEN
        
        energy(this_box)%ewald_reciprocal = E_reciprocal_move
        energy(this_box)%ewald_self = energy(this_box)%ewald_self + E_self_move
        
     END IF

     IF ( int_vdw_sum_style(this_box) == vdw_cut_tail) THEN
        energy(this_box)%lrc = e_lrc
     END IF

     ! obtain the original position of the deleted molecule so that the linked list
     ! can be updated

     do is = 1, nspecies
     k = 0
     CALL Get_Position_Molecule(this_box,is,im(is),position)

     IF (position < SUM(nmols(is,:))) THEN
        DO k = position + 1, SUM(nmols(is,:))
           locate(k-1,is) = locate(k,is)
        END DO
     END IF
     
     ! move the deleted molecule to the end of alive(is) molecules

     locate(nmols(is,this_box),is) = alive(is)
        
     molecule_list(alive(is),is)%live = .FALSE.
     atom_list(:,alive(is),is)%exist = .FALSE.

     nmols(is,this_box) = nmols(is,this_box) - 1
     enddo
     nsuccess(is,this_box)%deletion = nsuccess(is,this_box)%deletion + 1

!     CALL System_Energy_Check(1,mcstep,randno)
  ELSE

     is = 1
     IF ( (int_charge_sum_style(this_box) == charge_ewald) .AND. (has_charge(is)) ) THEN
        ! restore cos_sum and sin_sum. Note that these were changed when difference in
        ! reciprocal energies was computed

        !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)        
        cos_sum(:,this_box) = cos_sum_old(:,this_box)
        sin_sum(:,this_box) = sin_sum_old(:,this_box)
        !$OMP END PARALLEL WORKSHARE

     END IF

     IF ( int_vdw_sum_style(this_box) == vdw_cut_tail ) THEN
        ! restore the total number of bead types
        nint_beads(:,this_box) = nbeads_out(:)
     END IF

  END IF

!  IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)

  pacc = -pacc
  IF( pacc < 0.0_DP) THEN
    pacc = DEXP(pacc)
  ELSE
    pacc = 1.0_DP
  END IF

END SUBROUTINE Deletion

     

