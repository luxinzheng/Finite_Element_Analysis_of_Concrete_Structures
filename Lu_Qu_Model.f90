附录 A8 陆新征－曲哲滞回模型子程序

!=================================================================
!	陆新征－曲哲滞回模型
!	程序开发人：陆新征，陆新征
!	清华大学土木工程系
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
	gama= props(5)      ! 捏拢点荷载
	esoft=props(6)      ! 软化比例
	alpha= props(7)     ! 极限荷载和屈服荷载的比例
 	beta = props(8)     ! 正向与反向屈服强度比
	a_k= props(9)       ! 卸载刚度系数
	Omega= props(10)    ! 裂缝闭合位置, 1 很后 0 很早, Omega小于零为平行四边形模型
     
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

	if(a_k<0.) a_k=0.           ! 防止出现不当的计算结果 
	if(eta<=0.) eta=1.d-6;      ! 防止出现不当的计算结果
	if(esoft>=0.) esoft=-1.d-6  ! 防止出现不当的计算结果
      
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
		if(s0>0.) E_unload=E0*(abs(emax/(sy0/E0)))**(-a_k)  ! 计算卸载刚度
		if(s0<0.) E_unload=E0*(abs(emin/(sy0/E0)))**(-a_k)  ! 计算卸载刚度 
		if(E_unload<0.1*E0) E_unload=0.1*E0
	else
		E_unload=E0
	end if

	if(s0*de<0.) then ! 卸载行为
		s = s0 + E_unload * de
		Et = E_unload
		if(s*s0<0.) then ! 荷载出现反向, 开始加载或者再加载
			de=de-s0/E_unload
			s0=1.D-6*sy0*sign(1.d0,s) ! 给应力赋予一很小值
			Et=E0
		end if
	end if

		
      if ( de .ge. 0.0 .and. s0>=0.) then ! 正向加载
          sy = (1.0 - dt) * sy0
          !loading envelope
		! 强化段
		if(e+de>sy/E0) then
			evs = max( sy + ( e + de - sy/E0) * eta * E0, 0.)
			evE = eta * E0
	       if (s .ge. evs) then
		      s = evs;  Et = evE
		   end if
		end if
		! 软化段
		epeak=sy/E0+(alpha-1.)*sy/E0/eta
		if(e+0.5*de>epeak) then
			evs=max(sy*alpha+esoft*E0*(e+de-epeak),0.)
			evE=esoft*E0
		    if (s .ge. evs) then
			   s = evs;   Et = evE
		  end if
		end if

          !reloading envelope
          smax = max(sy, sy + (emax - sy/E0) * eta * E0)  ! 更新smax
		if(emax>epeak)	then ! 如果荷载进入软化段
			smax=max(sy*alpha+esoft*E0*(emax-epeak),0.)
		end if	
          sres = 0.02 * smax !0.2 * smax                              ! 得到荷载误差判别准则
          eres = ert - (srt - sres) / E_unload                  ! 得到变形判别准则

		if(Omega>=0) then
			x=emax-smax/E0  ! 最大荷载卸载降低到零时对应的变形
			e_slip=gama*emax+(1.-gama)*x ! 滑移捏拢终点对应的变形
			s_slip=smax*gama ! 滑移捏拢终点对应的荷载
			e_close=e_slip*Omega ! 裂缝闭合点
			s_close=(e_close-eres)/(e_slip-eres) * (s_slip-sres) + sres
		else
			e_slip=eres+sy*gama/E_unload
			s_slip=smax*gama
			e_close=e_slip; s_close=s_slip
		end if

          if (eres .le. emax - smax / E0) then    ! 曾经发生过反向加载
			if(e+0.5*de<e_close)  then  ! 变形小于卸载零点(裂缝闭合点), 此时发生捏拢
				srel = (e+de-eres)/(e_slip-eres) * (s_slip-sres) + sres
				Et1=(s_slip-sres)/(e_slip-eres)
			else                  ! 变形大于裂缝闭合点，指向历史最大点
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
		! 强化段
		if(e+de<-sy/E0) then
			evs =  min(-sy + ( e + de + sy/E0) * eta * E0,0.)
			evE = eta * E0
			if (s .le. evs) then
				s = evs; Et = evE
			end if
		end if
		! 软化段
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
 		if(emin<epeak)	then ! 如果荷载进入软化段
			smin=min(-sy*alpha+esoft*E0*(emin-epeak),0.)
		end if	
          sres = 0.02 * smin ! 0.2 * smin
          eres = erc - (src - sres) /  E_unload

		if(Omega>=0) then
			x=emin-smin/E0  ! 最大荷载卸载降低到零时对应的变形
			e_slip=gama*emin+(1.-gama)*x ! 滑移捏拢终点对应的变形
			s_slip=smin*gama ! 滑移捏拢终点对应的荷载
			e_close=e_slip*Omega ! 裂缝闭合点
			s_close=(e_close-eres)/(e_slip-eres) * (s_slip-sres) + sres
		else
			e_slip=eres-sy*gama/E_unload
			s_slip=smin*gama
			e_close=e_slip;	s_close=s_slip
		end if

          if (eres .ge. emin - smin / E0) then    ! 曾经发生过反向加载
			if(e+0.5*de>e_close) then ! 变形小于卸载零点(裂缝闭合点), 此时发生捏拢
				srel = (e+de-eres)/(e_slip-eres) * (s_slip-sres) + sres
				Et1=(s_slip-sres)/(e_slip-eres)
			else                  ! 变形大于裂缝闭合点，指向历史最大点
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