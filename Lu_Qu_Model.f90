��¼ A8 ½�����������ͻ�ģ���ӳ���

!=================================================================
!	½�����������ͻ�ģ��
!	���򿪷��ˣ�½������½����
!	�廪��ѧ��ľ����ϵ
!	 2009.10
!=================================================================

	subroutine Lu_Qu_Cycle(props,s,e,de,Et,statev,spd)
c
	IMPLICIT REAL *8 (A-H, O-Z)
	real*8 mu
	integer kon
c
      dimension props(10), statev(11)

      E0  = props(1)      ! Initial tangent stiffnss
	sy0 = props(2)      ! Initial yield stress
	eta = props(3)      ! Strain hardening ratio
	mu  = props(4)      ! maximum ductility
	gama= props(5)      ! ��£�����
	esoft=props(6)      ! ������
	alpha= props(7)     ! ���޺��غ��������صı���
 	beta = props(8)     ! �����뷴������ǿ�ȱ�
	a_k= props(9)       ! ж�ظն�ϵ��
	Omega= props(10)    ! �ѷ�պ�λ��, 1 �ܺ� 0 ����, OmegaС����Ϊƽ���ı���ģ��
     
      emax  = statev(1)   !maximum strain
      emin  = statev(2)   !minimum strain
      ert   = statev(3)   !stain at load reversal toward tension
      srt   = statev(4)   !stress at load reversal toward tension
      erc   = statev(5)   !stain at load reversal toward compression
      src   = statev(6)   !stress at load reversal toward compression
      kon   = nint(statev(7)) !
      Ehc   = statev(8)   !effective cummulative hysteresis energy
      Eh1   = statev(9)   !hysteresis energy in a half cycle
      dt    = statev(10)  !damage index for tension
      dc    = statev(11)  !damage index for compression
      eu    = mu * sy0/E0 !characteristic ultimate strain

	if(a_k<0.) a_k=0.           ! ��ֹ���ֲ����ļ����� 
	if(eta<=0.) eta=1.d-6;      ! ��ֹ���ֲ����ļ�����
	if(esoft>=0.) esoft=-1.d-6  ! ��ֹ���ֲ����ļ�����
      
      if (kon.eq.0) then
        emax =  sy0/E0
        emin = -beta*sy0/E0
        if (de.ge.0.0) then
            kon = 1
        else
            kon = 2
        end if
      else if ((kon.eq.1).and.(de.lt.0.0)) then !Load reversal
            kon = 2
            if (s.gt.0.0) then
                erc = e
                src = s
            end if
            Ehc = Ehc + Eh1 * (erc / eu ) ** 2.0
            Eh1 = 0.0 !a new half cycle is to begin
            if (e.gt.emax) emax = e
      else if ((kon.eq.2).and.(de.gt.0.0)) then !Load reversal
            kon = 1
            if (s.lt.0.0) then
                ert = e
                srt = s
            endif
            Ehc = Ehc + Eh1 * (ert / eu ) ** 2.0
            Eh1 = 0.0 !a new half cycle is to begin
            if (e.lt.emin) emin = e
      end if
c
	s0=s
	s = s + E0 * de
	Et = E0
	
	if(a_k>0.) then
		if(s0>0.) E_unload=E0*(abs(emax/(sy0/E0)))**(-a_k)  ! ����ж�ظն�
		if(s0<0.) E_unload=E0*(abs(emin/(sy0/E0)))**(-a_k)  ! ����ж�ظն� 
		if(E_unload<0.1*E0) E_unload=0.1*E0
	else
		E_unload=E0
	end if

	if(s0*de<0.) then ! ж����Ϊ
		s = s0 + E_unload * de
		Et = E_unload
		if(s*s0<0.) then ! ���س��ַ���, ��ʼ���ػ����ټ���
			de=de-s0/E_unload
			s0=1.D-6*sy0*sign(1.d0,s) ! ��Ӧ������һ��Сֵ
			Et=E0
		end if
	end if

		
      if ( de .ge. 0.0 .and. s0>=0.) then ! �������
          sy = (1.0 - dt) * sy0
          !loading envelope
		! ǿ����
		if(e+de>sy/E0) then
			evs = max( sy + ( e + de - sy/E0) * eta * E0, 0.)
			evE = eta * E0
	       if (s .ge. evs) then
		      s = evs;  Et = evE
		   end if
		end if
		! ����
		epeak=sy/E0+(alpha-1.)*sy/E0/eta
		if(e+0.5*de>epeak) then
			evs=max(sy*alpha+esoft*E0*(e+de-epeak),0.)
			evE=esoft*E0
		    if (s .ge. evs) then
			   s = evs;   Et = evE
		  end if
		end if

          !reloading envelope
          smax = max(sy, sy + (emax - sy/E0) * eta * E0)  ! ����smax
		if(emax>epeak)	then ! ������ؽ�������
			smax=max(sy*alpha+esoft*E0*(emax-epeak),0.)
		end if	
          sres = 0.02 * smax !0.2 * smax                              ! �õ���������б�׼��
          eres = ert - (srt - sres) / E_unload                  ! �õ������б�׼��

		if(Omega>=0) then
			x=emax-smax/E0  ! ������ж�ؽ��͵���ʱ��Ӧ�ı���
			e_slip=gama*emax+(1.-gama)*x ! ������£�յ��Ӧ�ı���
			s_slip=smax*gama ! ������£�յ��Ӧ�ĺ���
			e_close=e_slip*Omega ! �ѷ�պϵ�
			s_close=(e_close-eres)/(e_slip-eres) * (s_slip-sres) + sres
		else
			e_slip=eres+sy*gama/E_unload
			s_slip=smax*gama
			e_close=e_slip; s_close=s_slip
		end if

          if (eres .le. emax - smax / E0) then    ! �����������������
			if(e+0.5*de<e_close)  then  ! ����С��ж�����(�ѷ�պϵ�), ��ʱ������£
				srel = (e+de-eres)/(e_slip-eres) * (s_slip-sres) + sres
				Et1=(s_slip-sres)/(e_slip-eres)
			else                  ! ���δ����ѷ�պϵ㣬ָ����ʷ����
				srel=(e+de-e_close)/(emax-e_close)*(smax-s_close)+s_close
				Et1 = (smax - s_close) / (emax - e_close)
			end if
            if (s .gt. srel) then
               s = max( srel, 0.)
               Et = Et1
            end if
          end if

      elseif ( de .lt. 0.0 .and. s0<0. ) then
          sy = (1.0 - dc) * sy0 *beta
          !loading envelope
		! ǿ����
		if(e+de<-sy/E0) then
			evs =  min(-sy + ( e + de + sy/E0) * eta * E0,0.)
			evE = eta * E0
			if (s .le. evs) then
				s = evs; Et = evE
			end if
		end if
		! ����
		epeak=-sy/E0-(alpha-1.)*sy/E0/eta
		if(e+0.5*de<epeak) then
			evs=min(-sy*alpha+esoft*E0*(e+de-epeak),0.)
			evE=esoft*E0
		    if (s .le. evs) then
			    s = evs; Et = evE
		  end if
		end if

          !reloading envelope
          smin = min(-sy, -sy + (emin + sy/E0) * eta * E0)
 		if(emin<epeak)	then ! ������ؽ�������
			smin=min(-sy*alpha+esoft*E0*(emin-epeak),0.)
		end if	
          sres = 0.02 * smin ! 0.2 * smin
          eres = erc - (src - sres) /  E_unload

		if(Omega>=0) then
			x=emin-smin/E0  ! ������ж�ؽ��͵���ʱ��Ӧ�ı���
			e_slip=gama*emin+(1.-gama)*x ! ������£�յ��Ӧ�ı���
			s_slip=smin*gama ! ������£�յ��Ӧ�ĺ���
			e_close=e_slip*Omega ! �ѷ�պϵ�
			s_close=(e_close-eres)/(e_slip-eres) * (s_slip-sres) + sres
		else
			e_slip=eres-sy*gama/E_unload
			s_slip=smin*gama
			e_close=e_slip;	s_close=s_slip
		end if

          if (eres .ge. emin - smin / E0) then    ! �����������������
			if(e+0.5*de>e_close) then ! ����С��ж�����(�ѷ�պϵ�), ��ʱ������£
				srel = (e+de-eres)/(e_slip-eres) * (s_slip-sres) + sres
				Et1=(s_slip-sres)/(e_slip-eres)
			else                  ! ���δ����ѷ�պϵ㣬ָ����ʷ����
				srel=(e+de-e_close)/(emin-e_close)*(smin-s_close)+s_close
				Et1 = (smin - s_close) / (emin - e_close)
			end if
            if (s .lt. srel) then
                s = min (srel, 0.)
                Et = Et1
            end if
          end if
      end if


      if (Et.ne.E0 .and. Et.ne. E_unload) then 
            spd = spd + s * de
            Eh1 = Eh1 + s * de
            if ( s .ge. 0.0 ) then
                dc = min(Ehc /(3.0 * beta* sy0 * eu), 0.7)
            else
                dt = min(Ehc /(3.0 * sy0 * eu), 0.7)
            end if
      end if
c
	statev(1)   = emax
	statev(2)   = emin
	statev(3)   = ert
	statev(4)   = srt
	statev(5)   = erc
	statev(6)   = src
	statev(7)   = kon
	statev(8)   = Ehc
	statev(9)   = Eh1
	statev(10)  = dt
	statev(11)  = dc
	return
	end subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	