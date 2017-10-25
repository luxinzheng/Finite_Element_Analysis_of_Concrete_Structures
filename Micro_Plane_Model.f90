附录 A9 第四代微平面模型子程序


!=================================================================
!	本程序基于Bazant and Caner (2000)提出的第四代微平面模型
!	程序开发人：陆新征，黄羽立
!	清华大学土木工程系
!	 2007.3
!=================================================================

!*************
	module NumKind
!*************
	implicit none
		integer (kind(1)), parameter	:: ikind = kind(1),
	1		 rkind = kind(0.D0), lkind = kind(.true.)
	end module Numkind

!*********************
	module MicroplaneParam
!*********************
	use NumKind
	implicit none
	
	! Constants
	! Tolerance
	real (rkind),      parameter	:: Toler = 1.E-6
	! Number of microplane
	integer (ikind), parameter :: nMicroplane = 28
 	! Table of normals
	real (rkind), parameter :: N(nMicroplane,3) = reshape( 
	1(/ 0.577350258827209,  0.577350258827209,  0.577350258827209,
	1   0.577350258827209,  0.935113131999969,  0.935113131999969,
	1   0.935113131999969,  0.935113131999969,  0.250562787055969,
	1   0.250562787055969,  0.250562787055969,  0.250562787055969,  
	1   0.250562787055969,  0.250562787055969,  0.250562787055969,  
	1   0.250562787055969,  0.186156719923019,  0.186156719923019,  
	1   0.186156719923019,  0.186156719923019,  0.694746613502502,  
	1   0.694746613502502,  0.694746613502502,  0.694746613502502,  
	1   0.694746613502502,  0.694746613502502,  0.694746613502502,  
	1   0.694746613502502,  0.577350258827209,  0.577350258827209,  
	1  -0.577350258827209, -0.577350258827209,  0.250562787055969,  
	1   0.250562787055969, -0.250562787055969, -0.250562787055969,  
	1   0.935113131999969,  0.935113131999969, -0.935113131999969,  
	1  -0.935113131999969,  0.250562787055969,  0.250562787055969,
	1  -0.250562787055969, -0.250562787055969,  0.694746613502502,  
	1   0.694746613502502, -0.694746613502502, -0.694746613502502,  
	1   0.186156719923019,  0.186156719923019, -0.186156719923019,  
	1  -0.186156719923019,  0.694746613502502,  0.694746613502502,  
	1  -0.694746613502502, -0.694746613502502,  0.577350258827209,  
	1  -0.577350258827209,  0.577350258827209, -0.577350258827209,  
	1   0.250562787055969, -0.250562787055969,  0.250562787055969,  
	1  -0.250562787055969,  0.250562787055969, -0.250562787055969,  
	1   0.250562787055969, -0.250562787055969,  0.935113131999969,  
	1  -0.935113131999969,  0.935113131999969, -0.935113131999969,  
	1   0.694746613502502, -0.694746613502502,  0.694746613502502,  
	1  -0.694746613502502,  0.694746613502502, -0.694746613502502,  
	1   0.694746613502502, -0.694746613502502,  0.186156719923019,
	1  -0.186156719923019,  0.186156719923019, -0.186156719923019/),  
	1    (/nMicroplane,3/) )

	! Table of weights
	real (rkind), parameter :: w(nMicroplane) = 
	1 (/0.016071427613, 0.016071427613, 0.016071427613,  
	1   0.016071427613, 0.020474473014, 0.020474473014,  
	1   0.020474473014, 0.020474473014, 0.020474473014,  
	1   0.020474473014, 0.020474473014, 0.020474473014,  
	1   0.020474473014, 0.020474473014, 0.020474473014,  
	1   0.020474473014, 0.015835050493, 0.015835050493,  
	1   0.015835050493, 0.015835050493, 0.015835050493,  
	1   0.015835050493, 0.015835050493, 0.015835050493,  
	1   0.015835050493, 0.015835050493, 0.015835050493,  
	1   0.015835050493/)                               

 	! Fixed parameters
	real (rkind), parameter :: c1  = 0.62E0, 
	1 c2  = 2.76E0, 
	1 c3  = 4.E0,   
	1 c4  = 7.E1,   
	1 c5  = 2.5E0,  
	1 c6  = 1.3E0,  
	1 c7  = 50.E0,  
	1 c8  = 8.E0,   
	1 c9  = 1.3E0,  
	1 c10 = 0.73E0, 
	1 c11 = 0.2E0,  
	1 c12 = 7.E3,   
	1 c13 = 0.2E0,  
	1 c14 = 0.2E0,  
	1 c15 = 0.02E0, 
	1 c16 = 0.01E0, 
	1 c17 = 0.4E0,  
	1 c18 = 0.12E0

	real (rkind) :: k1, k2, k3, k4
	! Kronecker Delta
	real (rkind), parameter :: Kronecker(6) = (/1.E0, 1.E0, 1.E0, 
	1	0.E0, 0.E0, 0.E0/)

	real (rkind) :: E, nu, EV, ED, ET
	real (rkind) :: De(6,6)
	! Nij, Mij, Lij
	real (rkind) :: Nij(nMicroplane,6), Mij(nMicroplane,6), 
	1	Lij(nMicroplane,6)
	! Flag of initilization
	logical (lkind) :: bInitialized = .false.

	contains

	!===================================
	subroutine Initialize(props, nprops)
	!===================================

		integer (ikind), intent (in) :: nprops
		real    (rkind), intent (in) :: props(nprops)

		! Vector to tensor
		integer (ikind), parameter :: V2T(6,2) = reshape
	1		((/1,2,3,1,1,2,1,2,3,2,3,3/), (/6,2/))

		integer (ikind) :: i, j, t
		! Vector M, L
		real (rkind) :: M(3), L(3)
		real (rkind) :: len
		real (rkind) :: K11, K12, G

		!********************************************
!		!	props(1) - compressive strength
		!	props(1) - Elastic Modulus
		!	props(2) - Poisson's ratio
		!	props(3) - k1
		!	props(4) - k2
		!	props(5) - k3
		!	props(6) - k4
		!********************************************
		E = props(1)
		nu = props(2)
		k1 = props(3)
		k2 = props(4)
		k3 = props(5)
		k4 = props(6)

		if (nu > 0.5E0 .or. nu < 0.E0) then
			nu = 0.12E0
		end if

		EV = E / (1.E0 - 2.E0 * nu)
		ED = E / (1.E0 + nu)
		ET = ED

		do i = 1, nMicroplane

			len = sqrt(N(i,1) ** 2 + N(i,2) ** 2)
			if (abs(len) > Toler) then
				M = (/N(i,2) / len, -N(i,1) / len, 0.D0/)
			else
				M = (/1.E0, 0.E0, 0.E0/)
			end if

			L = (/ M(2) * N(i,3) - M(3) * N(i,2), 
	1			   M(3) * N(i,1) - M(1) * N(i,3), 
	1			   M(1) * N(i,2) - M(2) * N(i,1) /)

			do j = 1, 6
				Nij(i,j) = N(i,V2T(j,1)) * N(i,V2T(j,2))
				Mij(i,j) = 0.5E0 * ( M(V2T(j,1)) * N(i,V2T(j,2)) + 
	1				M(V2T(j,2)) * N(i,V2T(j,1)) )
				Lij(i,j) = 0.5E0 * ( L(V2T(j,1)) * N(i,V2T(j,2)) + 
	1				L(V2T(j,2)) * N(i,V2T(j,1)) )
			end do

		end do

		G = ED / 2.E0
		K12 = E * nu / (1.E0 + nu) / (1.E0 - 2.E0 * nu)
		K11 = K12 + ED

		! elastic stiffness matrix
		De(:,1) = (/ K11,  K12,  K12, 0.D0, 0.D0, 0.D0/)
		De(:,2) = (/ K12,  K11,  K12, 0.D0, 0.D0, 0.D0/)
		De(:,3) = (/ K12,  K12,  K11, 0.D0, 0.D0, 0.D0/)
		De(:,4) = (/0.D0, 0.D0, 0.D0,   ED, 0.D0, 0.D0/)
		De(:,5) = (/0.D0, 0.D0, 0.D0, 0.D0,   ED, 0.D0/)
		De(:,6) = (/0.D0, 0.D0, 0.D0, 0.D0, 0.D0,   ED/)

		bInitialized = .true.

	end subroutine Initialize

	end module MicroplaneParam

!********************************************************************
	subroutine CalStress(stress, statev, strain0, dstrain, ntens, nstatv)
!********************************************************************
	use NumKind
	use MicroplaneParam
	implicit none

	!****************************************************************
	!	statev(                                        1) - SigVpre
	!	statev(                  2 :     nMicroplane + 1) - SigNpre
	!	statev(    nMicroplane + 2 : 2 * nMicroplane + 1) - SigMpre
	!	statev(2 * nMicroplane + 2 : 3 * nMicroplane + 1) - SigLpre
	!****************************************************************
	integer (ikind), intent (in    ) :: ntens, nstatv
	real    (rkind), intent (in out) :: stress(ntens), statev(nstatv)
	real    (rkind), intent (in    ) :: strain0(ntens), dstrain(ntens)

	integer (ikind) :: i, j
	real (rkind) :: CV, CD, CT
	real (rkind) :: EpsV0, EpsD0, EpsT0, EpsN0
	real (rkind) :: EpsN, EpsV, EpsD, EpsM, EpsL
	real (rkind) :: dEpsN, dEpsV, dEpsD, dEpsM, dEpsL
	real (rkind) :: SigVpre, SigDpre, SigMpre, SigLpre, SigTpre,sigNpre
	real (rkind) :: SigVe, SigDe, SigMe, SigLe, SigTe
	real (rkind) :: SigNb, SigVbneg, SigVbpos, SigDbneg, SigDbpos, SigTb
	real (rkind) :: SigVstar, SigV, SigD, SigN, SigM, SigL
	real (rkind) :: SigN0, rateSigT, SumSigN
	real (rkind) :: strain1(ntens)

	real (rkind) :: EpsDArray(nMicroplane)

	strain1 = strain0 + dstrain

	!Volumetric microstrain
	EpsV = (strain1(1) + strain1(2) + strain1(3)) / 3.E0
	dEpsV = (dstrain(1) + dstrain(2) + dstrain(3)) / 3.E0
	SigVpre = statev(1)

	!Volumetric microstress
	EpsV0 = EpsV - dEpsV
	if (EpsV0 <= 0.E0 .and. SigVpre <= 0.E0 .and. dEpsV > 0.E0) then
		CV = EV * (c15 / (c15 - EpsV0) + SigVpre * EpsV0 / (c15 * c16 * EV))
	elseif (EpsV0 > 0.E0 .and. SigVpre > 0.E0 .and. dEpsV < 0.E0) then
		CV = min(SigVpre / EpsV0, EV)
	else
		CV = max(EV, E * k3 / k4 * exp(-EpsV0 / k1 / k4))
	end if

	SigVe = SigVpre + CV * dEpsV
	SigVbneg = -E * k1 * k3 * exp(-EpsV / k1 / k4)
	SigVbpos = EV * k1 * c13 / (1.E0 + c14 / k1 * 
	1	max(EpsV - k1 * c13, 0.E0)) ** 2
	SigVstar = min(max(SigVe, SigVbneg), SigVbpos)

	!clear sumation of SigN
	SumSigN = 0.E0

	do i = 1, nMicroplane
		! macrostrain to microstrain

		! normal microstrain
		EpsN = 0.E0
		dEpsN = 0.E0
		do j = 1, 6
			if(j<4) then
			EpsN = EpsN + Nij(i,j) * strain1(j)
			dEpsN = dEpsN + Nij(i,j) * dstrain(j)
			else 
			EpsN = EpsN + Nij(i,j) * strain1(j)* real(2)
			dEpsN = dEpsN + Nij(i,j) * dstrain(j)*real(2)
			end if
		end do

		EpsN0=EpsN-dEpsN
		! diviatoric microstrain
		EpsD = EpsN - EpsV
		dEpsD = dEpsN - dEpsV
		! restore diviatoric microstress
		! SigDpre = SigNpre - SigVpre
		sigNpre = statev(i + 1)
		SigDpre = statev(i + 1) - SigVpre

		EpsDArray(i) = EpsD

		! microstrain to microstress
		EpsD0 = EpsD - dEpsD
		if (     (EpsD0 > 0.E0 .and. SigDpre > 0.E0 .and. dEpsD < 0.E0) 
	1		.or. (EpsD0 < 0.E0 .and. SigDpre < 0.E0 .and. dEpsD > 0.E0) ) then
			CD = (1.E0 - c17) * ED + c17 * min(SigDpre / EpsD0, ED)
		else
			CD = ED
		end if

		! calculate elastic deviatoric microstress
		SigDe = SigDpre + CD * dEpsD
		! calculate boundary stress
		SigDbneg = -E * k1 * c8 / (1.E0 + (max(-EpsD - c8 
	1		* c9 * k1, 0.E0) / (k1 * c7)) ** 2)
		SigDbpos =  E * k1 * c5 / (1.E0 + (max( EpsD - c5 
	1		* c6 * k1, 0.E0) / (k1 * c18 * c7)) ** 2)
		! update deviatoric microstress
		SigD = min(max(SigDe, SigDbneg), SigDbpos)

		! normal microstress
		! calculate elastic deviatoric microstress
		SigN = SigVstar + SigD
		! calculate boundary stress
		SigNb = E * k1 * c1 * 
	1		exp( -max((EpsN - c1 * c2 * k1), 0.E0) / 
     1		(k1 * c3 + max(-c4 * SigVpre / EV, 0.E0)) )
		! update normal microstress
		SigN = min(SigN, SigNb)

		! crack-closing boundary
	IF (  EpsN0 * dEpsN< 0.e0 ) THEN
     	if (EpsN > 0.E0 .and. SigN < 0.E0) SigN = 0.E0
	ENDIF

		! sumation of normal microstress
		SumSigN = SumSigN + SigN * w(i)

		! store normal microstress
		statev(i + 1) = SigN

	end do

	SigV = min(SumSigN * 2.E0, SigVstar)

	statev(1) = SigV

	stress = 0.E0

	do i = 1, nMicroplane

		! restore SigN
		SigN =  statev(i + 1)

		! recalcualte Deviatoric
		SigD = SigN - SigV

		! macrostrain to microstrain
		EpsM = 0.E0
		dEpsM = 0.E0
		EpsL = 0.E0
		dEpsL = 0.E0
		do j = 1, 6
			if(j<4) then			
			EpsM = EpsM + Mij(i,j) * strain1(j)
			dEpsM = dEpsM + Mij(i,j) * dstrain(j)
			EpsL = EpsL + Lij(i,j) * strain1(j)
			dEpsL = dEpsL + Lij(i,j) * dstrain(j)
			else
			EpsM = EpsM + Mij(i,j) * strain1(j) * real(2)
			dEpsM = dEpsM + Mij(i,j) * dstrain(j)*real(2)
			EpsL = EpsL + Lij(i,j) * strain1(j)*real(2)
			dEpsL = dEpsL + Lij(i,j) * dstrain(j)*real(2)
			end if
		end do

		! restore microstress
		SigMpre = statev(i + nMicroplane + 1)
		SigLpre = statev(i + 2 * nMicroplane + 1)
		SigTpre = sqrt(SigMpre ** 2 + SigLpre ** 2)

		! check the loading criterion in shear
		EpsT0 = sqrt((EpsM - dEpsM) ** 2 + (EpsL - dEpsL) ** 2)
		if (sqrt(EpsM ** 2 + EpsL ** 2) - EpsT0 < 
	1		0.E0 .and. EpsT0 > Toler) then
			CT = (1.E0 - c17) * ET + c17 * min(SigTpre / EpsT0, ET)
!			write(7,*) 'T unloading in inc',kinc,'plane no.',i
		else
			CT = ET
		end if

		! shear microstress
		SigMe = SigMpre + CT * dEpsM
		SigLe = SigLpre + CT * dEpsL
		SigTe = sqrt(SigMe ** 2 + SigLe ** 2)
		SigN0 = ET * k1 * c11 / (1.E0 + c12 * max(EpsV, 0.E0))
		SigTb = ET * k1 * k2 * c10 * max(SigN0 - SigN, 0.E0) / 
	1		( ET * k1 * k2 + c10 * max(SigN0 - SigN, 0.E0) )

		if (abs(SigTe) > Toler) then
			rateSigT = min(1.E0, SigTb / SigTe)
		else
			rateSigT = 1.E0
		end if
		SigM = SigMe * rateSigT
		SigL = SigLe * rateSigT


		statev(i + nMicroplane + 1) = SigM
		statev(i + 2 * nMicroplane + 1) = SigL

		! microstress to macrostress
		do j = 1, 6
			stress(j) = stress(j) 
	1			+ ((Nij(i,j) - Kronecker(j) 
	1			/ 3.E0) * SigD + Mij(i,j) * SigM + Lij(i,j) * SigL)  * w(i)
		end do

	end do

	stress = stress * 6.E0
	stress(1:3) = stress(1:3) + SigV

	return
	end subroutine CalStress

