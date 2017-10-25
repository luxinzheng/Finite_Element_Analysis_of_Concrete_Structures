C	=====================================================
C	清华大学研究生课程《钢筋混凝土有限元》教学程序
C				等强双线性硬化本构模型				
C	利用MSC. Marc软件提供的 zero wkslp assoc yiel 等4个子程序
C	编制等强硬化双线性本构程序
C
C	主要参考文献：
C	1. MARC Volumn D user subroutine and Special Routines"
C	=====================================================
C	子程序说明：
C	zero:  得到等效应力
C	yiel:  得到当前屈服应力
C	wkslp: 得到强化模量
C	assoc: 得到流动方向

      subroutine assoc(stot,sinc,sc,t,ngens,ndi,nshear,n,nnn,layer)
c* * * * * *
c     provides flow rule for general plasticity option.
C		流动法则
C		变量说明
c     stoc         stress array !应力矩阵
c     sinc         flow direction array !流动向量
c     sc           hydrostatic stress !静水压力
c     t            not used
c     ngens        number of stress components ! 应力分量个数
c     ndi (3)       number of direct stress components ! 正应力分量个数
c     nshear (3)    number of shear stress components ! 剪应力分量个数
c     n            element number ! 单元编号
c     nnn          integration point number !积分点编号
c     layer        layer number ! 层编号，对于壳单元或者复合材料，可能有多个层
c* * * * * *
      implicit real*8 (a-h,o-z)   ! 变量说明
      dimension stot(*),sinc(*)  ! 变量说明
      do 26 i=1,ndi
26    sinc(i)=0.5*(3.*stot(i)-sc)  
      if(nshear.eq.0)go to 30
      f3=3.
      do 28 i=1,nshear
      j=ndi+i
28    sinc(j)=f3*stot(j)
30    continue
      return
      end

      subroutine wkslp(m,nn,kc,mats,slope,ebarp,eqrate,stryt,dt,ifirst)
      implicit real*8 (a-h,o-z)                                               dp
c* * * * * *
c     user subroutine to define work hardening slope.
c			得到硬化模量子程序
c     m            is the current user element number ! 单元编号
c     nn           is the integration point number ! 积分点编号
c     kc           is the layer number (=1 for non-layered elements) ! 层编号
c     mats         is the current material id ! 材料标号
c     slope        work hardening slope to be defined in this routine ! 返回硬化模量
c     ebarp        is the total equivalent plastic strain ! 等效塑性应变
c     eqrate       is the equivalent plastic strain rate ! 等效塑性应变率
c     stryt        current yield stress that optionally can be defined in this routine ! 返回
c                 屈服应力
c     dt           is the current total temperature ! 总温度
c     ifirst       flag distinguishing tenth cycle properties for
c                  ornl option
c
c    the internal element number mint can be obtained with
c      mint=ielint(m)
c* * * * * *
C	双线性弹塑性硬化
	SLOPE=20e3 !强化模量为Es=2GPa
	STRYT=210+ebarp*20e3 !强化后的屈服应力为210MPa+eps_eq_pl*Es
      return
      end

      real*8 function yiel(m,nn,kc,yld,ifirst,dt,eplas,erate,matz,
     *                        jprops)
c* * * * * *
C	得到屈服应力
c     find current value of yield stress
c     n            element index ! 单元编号
c     yld          initial yield stress !初始屈服应力
c     ifirst       flag for 10th cycle yield
c     dt           temperature ! 温度
c     eplas        equivalent plastic strain ! 等效塑性应变
c     erate        equivalent plastic strain rate !等效塑性应变率
c	jprops       table asscociated with the yield ! 定义屈服应力应变关系表（可选用）
c* * * * * *
      implicit real*8 (a-h,o-z)
      yiel=yld
c
c     user subroutine wkslp
c	调用做功强化子程序得到强化模量和屈服应力
         stryt=0.d0
         call wkslp(m,nn,kc,mats,slope,eplas,erate,stryt,dt,ifirst)
         yiel=stryt
      return
      end

      real*8 function zero(ndi,nshear,t,iort,ianiso,yrdir,yrshr,amm,a0)
c* * * * * *
c	得到等效应力
c     find equivalent stress for a stress tensor.
c     ndi        number of direct components of stress ! 正应力分量个数
c     nshear     number of shear  components of stress !剪应力分量个数
c     t(i)       component i of stress ! 应力向量
c     iort       indicating if curvilinear coordinates are used ! 是否使用曲线坐标系
c        iort = 0 - no curvilinear coordinates are used
c        iort = 1 - get mixed components in t
c        iort = 2 - assumes mixed compcomponents already in t
c     ianiso     flag indicating if anisotropy is used ! 是否为各向异性
c     yrdir      direct components for hill's anisotropic plasticity 
c     yrshr      shear  components for hill's anisotropic plasticity
c     amm        the metric if curvilinear coordinates are used
c     a0         the metric scale factor if curvilinear coordinates are used
c
c* * * * * *
      implicit real*8 (a-h,o-z)
      dimension t(*)
      dimension s(6),yrdir(3),yrshr(3),amm(3)
      zero=0.d0
c
c	等效应力采用von mises 应力
         if(ndi.eq.1) then
            zero=t(1)*t(1)*2.d0
	       else if(ndi.eq.2) then
               zero=2.d0*(t(1)*t(1)+t(2)*t(2)-t(1)*t(2))
                 else if(ndi.eq.3) then
                    zero=t(1)*t(1)+t(2)*t(2)+t(3)*t(3)
     *                 -t(1)*t(2)-t(2)*t(3)-t(3)*t(1)
                    zero=2.d0*zero
          endif
          if(nshear.eq.1) then
                zero=zero+6.d0*t(ndi+1)*t(ndi+1)
                else if(nshear.eq.2) then
                  zero=zero+6.d0*(t(ndi+1)*t(ndi+1)+t(ndi+2)*t(ndi+2))
                    else if(nshear.eq.3) then
                       zero=zero+6.d0*(t(ndi+1)*t(ndi+1)
     *                     +t(ndi+2)*t(ndi+2)+t(ndi+3)*t(ndi+3))
          endif
      if(zero.lt.0.d0) zero=0.d0
      zero=sqrt(0.5d0*zero)
      return
end


