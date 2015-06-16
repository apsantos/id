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

SUBROUTINE Insertion(this_box,mcstep,randno)

  !*************************************************************************
  ! 
  ! The subroutine inserts a molecule in the system. Various configurational
  ! schemes will be enabled as the code is verified.
  !
  ! Called by
  !
  !    gcmc_driver
  !
  ! Revision history
  !
  !   12/10/13 : Beta version
  !
  !*******************************************************************************

  USE Run_Variables
  USE Energy_Routines
  USE IO_Utilities
  USE Random_Generators
  USE Rotation_Routines
  USE Fragment_Growth

  IMPLICIT NONE

  INTEGER :: is, is_1, nmolecules_is, i, this_box, nmol_box, which_cell, i_type
  INTEGER :: nmols_sorbate, nsorbate, is_counter, rand_igas
  INTEGER, ALLOCATABLE :: sorbate_id(:), frag_order(:)
  INTEGER :: which_anchor,mcstep, ktothen, ifrag

  REAL(DP) :: dx, dy, dz, E_bond, E_angle, E_dihedral, E_intra_vdw, E_intra_qq
  REAL(DP) :: E_inter_vdw, E_inter_qq, E_improper
  REAL(DP) :: delta_e, E_reciprocal_move, E_self_move, E_lrc

  REAL(DP) :: factor, alpha_ratio, pick_species, randno
  REAL(DP) :: gnew, gold, dg  ! weighting factors for the new mol number and old mol number
  REAL(DP), ALLOCATABLE :: sorbate_x(:)
  REAL(DP) :: this_lambda

  LOGICAL :: accept, accept_or_reject 

  REAL(DP) :: checke, energy_old, energy_change
  LOGICAL  :: superbad

!FSL
  REAL(DP) :: tde, tintervdw, tinterqq, tb, ta, td, ti, ntemp, p_vdw, p_qq
  REAL(DP) :: tintravdw, tintraqq, tself, tpacc, suben, adden
  INTEGER :: tot_mols(nspecies), alive(nspecies), sm, am
  REAL(DP) :: P_forward(nspecies), nrg_ring_frag_tot(nspecies)
  LOGICAL :: cbmc_overlap(nspecies), inter_overlap(nspecies)
  LOGICAL :: intra_overlap(nspecies), p_overlap, rejectpair

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ntemp = 0.0_DP
  pacc = 0.0_DP
  paccbiased = 0.0_DP
  alpha_ratio = 1.0_DP
  inter_overlap(:) = .FALSE.
  cbmc_overlap(:) = .FALSE.
  P_forward(:) = 1.0_DP
  nrg_ring_frag_tot(:) = 0.0_DP
  tpacc = 0.0_DP
  p_vdw = 0.0_DP
  p_qq = 0.0_DP
  p_overlap = .FALSE.
  rejectpair = .FALSE.

  is = 0

  IF (int_sim_type .NE. sim_gemc) THEN
     this_box = INT(rranf() * nbr_boxes) + 1
     ! else for GEMC simulation the box will be specified as input
  END IF

  tot_trials(this_box) = tot_trials(this_box) + 1  
  energy_change = energy(1)%inter_vdw

  ! Choose a species to insert

!  is_1 = INT(rranf() * nspec_insert) + 1
!  is_counter = 0

!  DO is = 1, nspecies
!     IF(species_list(is)%int_species_type == int_sorbate) is_counter = is_counter + 1
!     IF(is_counter == is_1) EXIT
!  END DO

  nmolecules_is = 0

!!$  IF (int_sim_type == sim_gcmc) THEN
!!$     IF ( is == 1 ) t_collect = .TRUE.
!!$  END IF

  do is = 1, nspecies
  tot_mols(is) = SUM(nmols(is,:))

  IF (tot_mols(is) == (nmolecules(is))) RETURN
  enddo
  ! exit if we are attempting an insertion above maximum number

  ntrials(is,this_box)%insertion = ntrials(is,this_box)%insertion + 1
!  tot_trials(this_box) = tot_trials(this_box) + 1

 ! Assign a locate number for this molecule

  do is = 1, nspecies
  IF ( locate(tot_mols(is)+1,is) == 0 ) THEN
     locate(tot_mols(is)+1,is) = tot_mols(is) + 1
     ! otherwise we will use the locate number of a previously deleted molecule that has 
     ! been moved to the end of the array.
  END IF

  alive(is) = locate(tot_mols(is)+1,is)
  molecule_list(alive(is),is)%which_box = this_box
  this_lambda = 1.0_DP
  molecule_list(alive(is),is)%cfc_lambda = this_lambda
  molecule_list(alive(is),is)%molecule_type = int_normal
  if(is == 1) sm = alive(is)
  if(is == 2) am = alive(is)
  enddo

  ! Randomly insert the COM in the simulation box.
  
  do is = 1, nspecies
  IF (species_list(is)%fragment .AND. (species_list(is)%int_insert & 
       .NE. int_igas) ) THEN

     del_Flag = .FALSE.
     get_fragorder = .TRUE.
     ALLOCATE(frag_order(nfragments(is)))
     CALL Build_Molecule(alive(is),is,this_box,frag_order,this_lambda,which_anchor,P_forward(is),nrg_ring_frag_tot(is),cbmc_overlap(is))
     molecule_list(alive(is),is)%live = .TRUE.
     atom_list(:,alive(is),is)%exist = .TRUE.

     ktothen = 1

     IF (nfragments(is) /= 0 ) THEN

        ktothen = ktothen * kappa_ins

        IF (kappa_rot /= 0 ) THEN

           ktothen = ktothen * kappa_rot

        END IF

        IF (kappa_dih /= 0 ) THEN

           DO ifrag = 1, nfragments(is) - 1

              ktothen = ktothen * kappa_dih
              
           END DO

        END IF

     END IF

     P_forward(is) = P_forward(is) * REAL(ktothen, DP)

     DEALLOCATE(frag_order)

  ELSE
 
     atom_list(:,alive(is),is)%exist = .true.
     molecule_list(alive(is),is)%live = .true.
     IF(species_list(is)%int_insert == int_random) THEN
     
        ! COM of the species from the initial configuration is 

        molecule_list(alive(is),is)%xcom = species_list(is)%xcom
        molecule_list(alive(is),is)%ycom = species_list(is)%ycom
        molecule_list(alive(is),is)%zcom = species_list(is)%zcom

        
        atom_list(:,alive(is),is)%rxp = init_list(:,1,is)%rxp
        atom_list(:,alive(is),is)%ryp = init_list(:,1,is)%ryp
        atom_list(:,alive(is),is)%rzp = init_list(:,1,is)%rzp


     ELSE IF(species_list(is)%int_insert == int_igas) THEN

        rand_igas = (rranf() * n_igas(is)) + 1

        molecule_list(alive(is),is)%xcom = molecule_list_igas(rand_igas,is)%xcom
        molecule_list(alive(is),is)%ycom = molecule_list_igas(rand_igas,is)%ycom
        molecule_list(alive(is),is)%zcom = molecule_list_igas(rand_igas,is)%zcom

        atom_list(:,alive(is),is)%rxp = atom_list_igas(:,rand_igas,is)%rxp
        atom_list(:,alive(is),is)%ryp = atom_list_igas(:,rand_igas,is)%ryp
        atom_list(:,alive(is),is)%rzp = atom_list_igas(:,rand_igas,is)%rzp

     END IF   
     
     CALL Rotate_Molecule_Eulerian(alive(is),is)
     
     IF ( box_list(this_box)%int_box_shape == int_cubic ) THEN
        
        molecule_list(alive(is),is)%xcom = (rranf() - 0.5_DP) * box_list(this_box)%length(1,1)
        molecule_list(alive(is),is)%ycom = (rranf() - 0.5_DP) * box_list(this_box)%length(2,2)
        molecule_list(alive(is),is)%zcom = (rranf() - 0.5_DP) * box_list(this_box)%length(3,3)
        
     END IF
     
     ! Coordinates obtained are for the initial coordinates so translate it to the current position

     IF(species_list(is)%int_insert == int_random) THEN
     
        dx = molecule_list(alive(is),is)%xcom - species_list(is)%xcom
        dy = molecule_list(alive(is),is)%ycom - species_list(is)%ycom
        dz = molecule_list(alive(is),is)%zcom - species_list(is)%zcom

     ELSE

        dx = molecule_list(alive(is),is)%xcom - molecule_list_igas(rand_igas,is)%xcom
        dy = molecule_list(alive(is),is)%ycom - molecule_list_igas(rand_igas,is)%ycom
        dz = molecule_list(alive(is),is)%zcom - molecule_list_igas(rand_igas,is)%zcom

     END IF        
     
     atom_list(:,alive(is),is)%rxp = atom_list(:,alive(is),is)%rxp + dx
     atom_list(:,alive(is),is)%ryp = atom_list(:,alive(is),is)%ryp + dy
     atom_list(:,alive(is),is)%rzp = atom_list(:,alive(is),is)%rzp + dz

     ! Now compute the energy difference due to the insertion of this molecules
     ! make this molecule alive(is)

     molecule_list(alive(is),is)%live = .TRUE.
     molecule_list(alive(is),is)%which_box = this_box
  
     atom_list(:,alive(is),is)%exist = .TRUE.
     molecule_list(alive(is),is)%cfc_lambda = 1.0_DP

  END IF
  enddo

  ! compute the distance of atom farthest from COM
 
  delta_e = 0.0_DP
  E_inter_vdw = 0.0_DP
  E_inter_qq = 0.0_DP

  do is = 1, nspecies
  IF(.NOT. cbmc_overlap(is)) THEN

    CALL Fold_Molecule(alive(is),is,this_box)
    CALL Get_COM(alive(is),is)
    CALL Compute_Max_COM_Distance(alive(is),is)

    ! Intra molecule energy

     tde = 0.0_DP
     
     tintervdw = 0.0_DP
     tinterqq = 0.0_DP
     
     CALL Compute_Molecule_Nonbond_Inter_Energy(alive(is),is,tintervdw,tinterqq,inter_overlap(is))
     E_inter_vdw = E_inter_vdw + tintervdw
     E_inter_qq = E_inter_qq + tinterqq

  END IF
  enddo


       CALL Compute_Molecule_Pair_Interaction(sm,1,am,2,this_box, &
                  p_vdw,p_qq,p_overlap)

     E_inter_vdw = E_inter_vdw - p_vdw
     E_inter_qq = E_inter_qq - p_qq

  delta_e = E_inter_vdw + E_inter_qq

  do is = 1, nspecies
  IF (inter_overlap(is) .OR. cbmc_overlap(is)) THEN
     ! reject the move immediately
     rejectpair = .TRUE.
     EXIT
  END IF
  enddo

  if (rejectpair) then
  do is = 1, nspecies
     molecule_list(alive(is),is)%live = .FALSE.
     atom_list(:,alive(is),is)%exist = .FALSE.
  enddo
  RETURN
  endif

  E_bond = 0.0_DP
  E_angle = 0.0_DP
  E_dihedral = 0.0_DP
  E_improper = 0.0_DP

  do is = 1, nspecies
  tb = 0.0_DP
  ta = 0.0_DP
  td = 0.0_DP
  ti = 0.0_DP
  IF(species_list(is)%int_insert == int_random) THEN

     CALL Compute_Molecule_Bond_Energy(alive(is),is,tb)
     CALL Compute_Molecule_Angle_Energy(alive(is),is,ta)
     CALL Compute_Molecule_Dihedral_Energy(alive(is),is,td)
     CALL Compute_Molecule_improper_Energy(alive(is),is,ti)

  ELSE IF (species_list(is)%int_insert == int_igas) THEN

     tb = energy_igas(rand_igas,is)%bond
     ta = energy_igas(rand_igas,is)%angle
     td = energy_igas(rand_igas,is)%dihedral
     ti = energy_igas(rand_igas,is)%improper

  END IF
  E_bond = E_bond + tb
  E_angle = E_angle + ta
  E_dihedral = E_dihedral + td
  E_improper = E_improper + ti
  enddo

  delta_e = delta_e + E_bond + E_angle + E_dihedral + E_improper

  E_intra_vdw = 0.0_DP
  E_intra_qq = 0.0_DP
  do is = 1, nspecies
  tintravdw = 0.0_DP
  tintraqq = 0.0_DP
  CALL Compute_Molecule_Nonbond_Intra_Energy(alive(is),is,tintravdw,tintraqq,intra_overlap(is))
  E_intra_vdw = E_intra_vdw + tintravdw
  E_intra_qq = E_intra_qq + tintraqq
  enddo

!FSL Se intra_overlap = true, rejeitar inserir molecula
  do is = 1, nspecies
  IF (intra_overlap(is)) THEN
     rejectpair = .TRUE.
   endif
   enddo

   if (rejectpair) then
   do is = 1, nspecies
     molecule_list(alive(is),is)%live = .FALSE.
     atom_list(:,alive(is),is)%exist = .FALSE.
   enddo
   RETURN
   endif

  delta_e = delta_e + E_intra_vdw + E_intra_qq

  E_reciprocal_move = 0.0_DP
  E_self_move = 0.0_DP
  do is = 1, nspecies
  tself = 0.0_DP
  IF ( (int_charge_sum_style(this_box) == charge_ewald) .AND. (has_charge(is)) ) THEN

     CALL Ins_Pairs_Ewald_Reciprocal_Energy_Difference(alive(is),alive(is),is,this_box,int_insertion,E_reciprocal_move)
     CALL Compute_Ewald_Self_Energy_Difference(alive(is),is,this_box,int_insertion,tself)
     E_self_move = E_self_move + tself

  END IF
  enddo

  delta_e = delta_e + E_self_move + E_reciprocal_move - energy(this_box)%ewald_reciprocal

  IF (int_vdw_sum_style(this_box) == vdw_cut_tail) THEN
     ! increase number of integer beads

     nbeads_in = nint_beads(:,this_box)

     do is = 1, nspecies
     DO i = 1, natoms(is)
        i_type = nonbond_list(i,is)%atom_type_number
        nint_beads(i_type,this_box) = nint_beads(i_type,this_box) + 1
     END DO
     enddo

     CALL Compute_LR_correction(this_box,e_lrc)
     delta_e = delta_e + e_lrc - energy(this_box)%lrc

  END IF

  ! calculate weighting function and apply acceptance rule
  
    dg = 0.0_DP

    ntemp = 0.0_DP
    suben = 0.0_DP
    adden = 0.0_DP
!FSL energy    
  do is = 1, nspecies
  IF(species_list(is)%int_insert == int_igas) THEN 
     tpacc = beta(this_box) * (energy_igas(rand_igas,is)%total)
  ELSEIF (species_list(is)%fragment) THEN
     tpacc = beta(this_box) * (nrg_ring_frag_tot(is))
  ELSE
     tpacc = 0.0_DP
  END IF
  ntemp = ntemp + nmols(is,this_box) + 1
  if(is == 1) adden = adden + DLOG(P_forward(1)) + DLOG(P_forward(2))
  suben = suben + tpacc + beta(this_box)*E_angle*0.5_DP
  enddo
!FSL

  pacc = beta(this_box)*delta_e - suben + adden + DLOG(ntemp)

  pacc = pacc - DLOG(alpha_ratio) + DLOG(box_list(this_box)%volume)

  open(unit=58,file="ins_en.dat",action="write")
  write(58,*) "B1", E_inter_vdw, E_inter_qq 
  write(58,*) "B2", E_angle, E_dihedral
  write(58,*) "B3", E_intra_vdw, E_intra_qq
  write(58,*) "B4", E_self_move, E_reciprocal_move
  write(58,*) "B5", energy(this_box)%ewald_reciprocal, e_lrc, energy(this_box)%lrc
  write(58,*) "B6", pacc, delta_e, beta(this_box)*delta_e, &
     dlog(alpha_ratio), dlog(P_forward(1)), dlog(P_forward(2))

  write(58,*) "B7", dlog(dbpair), DLOG(species_list(1)%de_broglie(this_box)), &
     DLOG(species_list(2)%de_broglie(this_box))
  is = 1
  IF(lchempot) THEN
     pacc = pacc - species_list(1)%chem_potential * beta(this_box)
     pacc = pacc + 3.0_DP * DLOG(dbpair)
     pacc = pacc - DLOG (species_list(is)%zig_by_omega)
  ELSE
     pacc = pacc - DLOG(species_list(is)%fugacity * beta(this_box)) + &
                   DLOG(species_list(is)%zig_by_omega)
  END IF

  ! correct for ring biasing if any

  factor = pacc + dg  ! accept based on probability times weighting factor
 ! factor = pacc 
  
  accept = accept_or_reject(factor)
  
  factor = -factor
  IF( factor < 0.0_DP ) THEN
    factor = DEXP(factor)
  ELSE
    factor = 1.0_DP
  END IF

  paccbiased = factor

  IF (accept) THEN
     ! update the number of molecules
     do is = 1, nspecies
     nmols(is,this_box) = nmols(is,this_box) + 1
     enddo
     ! update the energies
     energy(this_box)%total = energy(this_box)%total + delta_e
     energy(this_box)%intra = energy(this_box)%intra + E_bond + E_angle + E_dihedral + E_improper
     energy(this_box)%bond = energy(this_box)%bond + E_bond
     energy(this_box)%angle = energy(this_box)%angle + E_angle
     energy(this_box)%dihedral = energy(this_box)%dihedral + E_dihedral
     energy(this_box)%improper = energy(this_box)%improper + E_improper
     energy(this_box)%intra_vdw = energy(this_box)%intra_vdw + E_intra_vdw
     energy(this_box)%intra_q = energy(this_box)%intra_q + E_intra_qq
     energy(this_box)%inter_vdw = energy(this_box)%inter_vdw + E_inter_vdw
     energy(this_box)%inter_q = energy(this_box)%inter_q + E_inter_qq

     IF ( int_charge_sum_style(this_box) == charge_ewald .AND. has_charge(is)) THEN

        energy(this_box)%ewald_reciprocal = E_reciprocal_move
        energy(this_box)%ewald_self = energy(this_box)%ewald_self + E_self_move
     END IF

     IF (int_vdw_sum_style(this_box) == vdw_cut_tail) THEN
        energy(this_box)%lrc = e_lrc
     END IF

     nsuccess(is,this_box)%insertion = nsuccess(is,this_box)%insertion + 1


  ELSE
 
     do is = 1, nspecies 
     molecule_list(alive(is),is)%live = .FALSE.
     atom_list(:,alive(is),is)%exist = .FALSE.
     molecule_list(alive(is),is)%molecule_type = int_none
     enddo

     is = 1 
     IF ( int_charge_sum_style(this_box) == charge_ewald .AND. (has_charge(is)) ) THEN
        ! restore cos_sum and sin_sum. Note that these were changed when difference in
        ! reciprocal energies was computed

        cos_sum(:,this_box) = cos_sum_old(:,this_box)
        sin_sum(:,this_box) = sin_sum_old(:,this_box)

     END IF

     IF ( int_vdw_sum_style(this_box) == vdw_cut_tail ) THEN

        nint_beads(:,this_box) = nbeads_in(:)

     END IF
     
  END IF

  ! Accept or reject the move.
  ! also, update various arrays that depend on the number of molecules
  pacc = -pacc
  IF( pacc < 0.0_DP ) THEN
    pacc = DEXP(pacc)
  ELSE
    pacc = 1.0_DP
  END IF
END SUBROUTINE Insertion
