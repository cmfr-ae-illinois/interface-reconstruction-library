!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2025 Ilia Kheirkhah <iliak2@illinois.edu>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

  
! In this example a quarter circle of radius 2.5
! is created using PLIC segments. These are then blended into
! an implicit surface by a Partition of Unity, which is then 
! used to calculate the forces on various edges. 
! Results are printed and verified.
program main
  use irl_fortran_interface
  use f_PUNeigh_RectCub_class
  use f_SeparatorVariant_class
  use f_PUSolve_RectCub_class
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  
  ! Declar Centroids
  real(DP), dimension(1:3) :: cen1,cen2,cen3,cen4,cen5,nor1,nor2,nor3,nor4,nor5
  real(DP), dimension(1:3) :: startPoint,endPoint,force,Marangoni
  real(DP) :: d1,d2,d3,d4,d5,stc,pressure
  real(DP),parameter :: delta = 5.0_DP
  type(SeparatorVariant_type) :: plane1,plane2,plane3,plane4,plane5
  ! First make a Neighborhood
  type(PUNeigh_RectCub_type) :: neighborhood
  ! Now the solver object
  type(PU_RectCub_type) :: solver

  ! Define Planar Separators
  ! Centroids
  write(*,'(A)') 'Centroids Making'
  cen1 = (/0.5_DP,2.4495097568_DP,0.0_DP /)
  cen2 = (/1.76776695297_DP,1.76776695297_DP,0.0_DP/)
  cen3 = (/1.26289662285_DP,2.15773797371_DP,0.0_DP/)
  cen4 = (/2.4495097568_DP,0.5_DP,0.0_DP/)
  cen5 = (/2.15773797371_DP,1.26289662285_DP,0.0_DP/)
  ! Normals
  write(*,'(A)') 'Normals Making'
  nor1 = (/0.196116135138_DP,0.980580675691_DP,0.0_DP/)
  nor2 = (/0.707106781187_DP,0.707106781187_DP,0.0_DP/)
  nor3 = (/0.514495755428_DP,0.857492925713_DP,0.0_DP/)
  nor4 = (/0.980580675691_DP,0.196116135138_DP,0.0_DP/)
  nor5 = (/0.857492925713_DP,0.514495755428_DP,0.0_DP/)
  ! Calculate Distances
  write(*,'(A)') 'Distances Making'
  d1 = nor1(1)*cen1(1) + nor1(2)*cen1(2) + nor1(3)*cen1(3)
  d2 = nor2(1)*cen2(1) + nor2(2)*cen2(2) + nor2(3)*cen2(3)
  d3 = nor3(1)*cen3(1) + nor3(2)*cen3(2) + nor3(3)*cen3(3)
  d4 = nor4(1)*cen4(1) + nor4(2)*cen4(2) + nor4(3)*cen4(3)
  d5 = nor5(1)*cen5(1) + nor5(2)*cen5(2) + nor5(3)*cen5(3)

  ! Make Separator Variants from normal and distance
  write(*,'(A)') 'Making Planes'
  call new(plane1)
  call setNumberOfPlanes(plane1,1)
  call setPlane(plane1, 0, nor1,d1)

  call new(plane2)
  call setNumberOfPlanes(plane2,1)
  call setPlane(plane2, 0, nor2,d2)

  call new(plane3)
  call setNumberOfPlanes(plane3,1)
  call setPlane(plane3, 0, nor3,d3)

  call new(plane4)
  call setNumberOfPlanes(plane4,1)
  call setPlane(plane4, 0, nor4,d4)

  ! Now, add Separators to Neighborhood
  write(*,'(A)') 'Making Neighborhood'
  call new(neighborhood)
  call addMember(neighborhood,cen1,1.0_DP,plane1)

  call new(plane1)
  call setNumberOfPlanes(plane1,1)
  call setPlane(plane1, 0, nor2,d2)
  call addMember(neighborhood,cen2,1.0_DP,plane1)

  ! call new(plane1)
  ! call setNumberOfPlanes(plane1,1)
  ! call setPlane(plane1, 0, nor3,d3)
  ! call addMember(neighborhood,cen3,1.0_DP,plane1)

  ! call new(plane1)
  ! call setNumberOfPlanes(plane1,1)
  ! call setPlane(plane1, 0, nor4,d4)
  ! call addMember(neighborhood,cen4,1.0_DP,plane1)

  ! call new(plane1)
  ! call setNumberOfPlanes(plane1,1)
  ! call setPlane(plane1, 0, nor5,d5)
  ! call addMember(neighborhood,cen5,1.0_DP,plane1)

  ! Now that everything is in the neighborhood, make the solver object and put the neighborhood in.
  write(*,'(A)') 'Making Solver'
  call new(solver)
  call setNeighborhood(solver,neighborhood)
  stc = 1.0_DP

  write(*,'(A)')
  write(*,'(A)') '============= RESULTS ============='
  ! The test ran is a quarter circle of radius 2.5. The intersection edge at x=0, between y=2 and y=3.
  ! The Points are calculated going form y=3 to y=2 to maintain counter-clockwise orientation.
  startPoint = (/0.0_DP,3.0_DP,0.0_DP/)
  endPoint = (/0.0_DP,2.0_DP,0.0_DP/)
  force = (/0.0_DP,0.0_DP,0.0_DP/)
  pressure = 0.0_DP
  Marangoni = (/0.0_DP,0.0_DP,0.0_DP/)
  call solveEdge(solver,stc,startPoint,endPoint,delta,pressure,Marangoni,force)
  write(*,'(A)') 'Quarter-Circle Force 1 Test' 
  write(*,'(A,3F10.5)') '> Expected: ', -0.961249020086_DP, 0.275681557931_DP, 0.0_DP
  write(*,'(A,3F10.5)') '> Computed: ', force(1),force(2),force(3)

  startPoint = (/1.0_DP,3.0_DP,0.0_DP/)
  endPoint = (/1.0_DP,2.0_DP,0.0_DP/)
  force = (/0.0_DP,0.0_DP,0.0_DP/)
  call solveEdge(solver,stc,startPoint,endPoint,delta,pressure,Marangoni,force)
  write(*,'(A)') 'Quarter-Circle Force 2 Test' 
  write(*,'(A,3F10.5)') '> Expected: ', -0.880916212182_DP, 0.473272254748_DP, 0.0_DP
  write(*,'(A,3F10.5)') '> Computed: ', force(1),force(2),force(3)

  startPoint = (/2.0_DP,2.0_DP,0.0_DP/)
  endPoint = (/1.0_DP,2.0_DP,0.0_DP/)
  force = (/0.0_DP,0.0_DP,0.0_DP/)
  call solveEdge(solver,stc,startPoint,endPoint,delta,pressure,Marangoni,force)
  write(*,'(A)') 'Quarter-Circle Force 3 Test' 
  write(*,'(A,3F10.5)') '> Expected: ', -0.739254963695_DP, 0.673425644487_DP, 0.0_DP
  write(*,'(A,3F10.5)') '> Computed: ', force(1),force(2),force(3)

  startPoint = (/2.0_DP,2.0_DP,0.0_DP/)
  endPoint = (/2.0_DP,1.0_DP,0.0_DP/)
  force = (/0.0_DP,0.0_DP,0.0_DP/)
  call solveEdge(solver,stc,startPoint,endPoint,delta,pressure,Marangoni,force)
  write(*,'(A)') 'Quarter-Circle Force 4 Test' 
  write(*,'(A,3F10.5)') '> Expected: ', -0.673425644487_DP, 0.739254963695_DP, 0.0_DP
  write(*,'(A,3F10.5)') '> Computed: ', force(1),force(2),force(3)

  startPoint = (/3.0_DP,1.0_DP,0.0_DP/)
  endPoint = (/2.0_DP,1.0_DP,0.0_DP/)
  force = (/0.0_DP,0.0_DP,0.0_DP/)
  call solveEdge(solver,stc,startPoint,endPoint,delta,pressure,Marangoni,force)
  write(*,'(A)') 'Quarter-Circle Force 5 Test' 
  write(*,'(A,3F10.5)') '> Expected: ', -0.473272254748_DP, 0.880916212182_DP, 0.0_DP
  write(*,'(A,3F10.5)') '> Computed: ', force(1),force(2),force(3)

  startPoint = (/3.0_DP,0.0_DP,0.0_DP/)
  endPoint = (/2.0_DP,0.0_DP,0.0_DP/)
  force = (/0.0_DP,0.0_DP,0.0_DP/)
  call solveEdge(solver,stc,startPoint,endPoint,delta,pressure,Marangoni,force)
  write(*,'(A)') 'Quarter-Circle Force 6 Test' 
  write(*,'(A,3F10.5)') '> Expected: ', -0.275681557931_DP, 0.961249020086_DP, 0.0_DP
  write(*,'(A,3F10.5)') '> Computed: ', force(1),force(2),force(3)


end program main
