!                                 Notice
!  This program was prepared by the University of California (University)
!  under Contract W-7405-ENG-36 with the U.S. Department of Energy (DOE).
!  All rights in the program are reserved by DOE on behalf of the Government
!  and the University pursuant to the contract. You are authorized to use
!  this program for Government purposes but it is not to be released or
!  distributed to the public.
!  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY,
!  NOR THE UNIVERSITY OF CALIFORNIA, NOR ANY OF THEIR EMPLOYEES,
!  MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY
!  OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS, OF
!  ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
!  THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
!
      module datamodule
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! module with all parameter and array specifications which are passed between main
! and subroutines
! created PKK 10/02/02
! arrays are allocatable, PKK 10/01/02
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         implicit none
         save
         integer ii,i,j,k,nx,ny,nz,inumbuild,ibuild,nz_data
         integer itermax, roofflag, upwindflag,intersectionflag,wakeflag,sidewall_flag ! MAN 6/29/2007
         integer, allocatable :: icellflag(:,:,:),ibldflag(:,:,:)
         integer, allocatable :: istart(:),iend(:)!e
         integer, allocatable :: jstart(:),jend(:)!e
         integer, allocatable :: kstart(:),kend(:)!e
         real, allocatable :: Lf(:),Lr(:)
         real, allocatable :: Weff(:),Leff(:)!e

         real, allocatable :: Lt(:)!e
!         real, allocatable :: Lfx(:),Lfy(:)!e
!         real, allocatable :: Lfx1(:),Lfy1(:) !NLB 02/10/04
         real, allocatable :: Rscale(:),Rcx(:) !NLB 10/11/04
!         real, allocatable :: Roofcx(:,:,:)  !NLB 10/11/04
         real, allocatable :: Wt(:)!e
         !MAN 08/10/2010 polygon buildings
         integer, allocatable :: bldnum(:),bldtype(:),bldgeometry(:)
         integer, allocatable :: bldstartidx(:),bldstopidx(:),numpolygons(:),bldroof(:)
         real, allocatable :: bldx(:),bldy(:),bldcx(:),bldcy(:),bldwall(:)
         !end MAN 08/10/2010
         real, allocatable :: uo_bl(:),vo_bl(:)	!erp 6/8/04
         real, allocatable :: uo(:,:,:),vo(:,:,:),wo(:,:,:)
         real, allocatable :: ufarwake(:,:),vfarwake(:,:)
         real, allocatable :: p1(:,:,:),p2(:,:,:),r(:,:,:)
         real, allocatable :: e(:,:,:),f(:,:,:),g(:,:,:)
         real, allocatable :: h(:,:,:),m(:,:,:),n(:,:,:)
         real, allocatable:: denom(:,:,:)
         real, allocatable :: u(:,:,:),v(:,:,:),w(:,:,:)
         real dx,dy,dz,Lx,Ly,Lz,omegarelax ! ,ddx MAN 7/7/2005 var dz conversion
         real A,B,eps,pi
         real alpha1,alpha2,eta
         real, allocatable :: Ht(:),Wti(:),Lti(:),aa(:),bb(:) !erp 1/31/2003
         real, allocatable :: xfo(:),yfo(:),zfo(:),gamma(:)   !erp 7/26/03
         real, allocatable :: ws_data(:),wd_data(:),z_data(:) !erp 6/5/2003
         real, allocatable :: u_data(:),v_data(:)			   !erp 6/5/2003
         real uin,zref,pp,zo ! MAN 8-19-2005 Updated input output file structures
         real logvel !NLB 04/20/04

         integer streetcanyonflag

!erp binary write
         real, allocatable :: x_write(:,:,:),y_write(:,:,:),z_write(:,:,:)
!erp
!multiple run variables
!         real start_time		 !decimal start time
!         real time_incr         !time increment
         integer, allocatable :: time(:) !time values in unix time
         integer num_time_steps,i_time

!variables added TMB 3/10/04
!time variable multi-sensor input
         real, allocatable :: site_xcoord(:), site_ycoord(:) !x,y coordinates of each sensor site
         real, allocatable :: u_prof(:,:), v_prof(:,:)						  !U and V vertical velocity profiles
         real, allocatable :: uoint(:), voint(:)							  !Interpolated U and V velocities
         real, allocatable :: wm(:,:,:),wms(:,:,:)								  !the "Weight" of each site on each grid point
         real, allocatable :: site_z_data(:,:,:)							  !z_data = height
         real, allocatable :: site_ws_data(:,:,:),site_wd_data(:,:,:)		  !ws_data = wind speed, wd_data = wind direction
         real, allocatable :: site_u_data(:,:,:),site_v_data(:,:,:)		  !u_data = u velocity, v_data = v velocity
         integer, allocatable :: site_nz_data(:,:),site_blayer_flag(:,:)  !number of points in the vertical direction and the boundary layer flag
         real, allocatable :: site_pp(:,:)                ! exp/zo for each site
         real, allocatable :: site_H(:,:),site_ac(:,:),site_rL(:,:)
! MAN 4/5/2007 unused variables: site_zcoord(:), t_site(:,:), dir_site(:,:), vel_site(:,:),site_ustar(:),site_lc(:)
         integer num_sites	!number of sensor sites
         integer num_vert_points !MAN 4/5/2007 number of data points to allocate in data points profiles
         integer theta_init  !the that is read each time step for blayers 1-4
!end variables added TMB 3/10/04
!erp file format flag 3/02/05
         integer format_flag,uofield_flag,uosensor_flag
!erp end file format flag
!erp 8/8/2005 adding uzma's veg model
         integer inumcanopy	!number of vegetative canopies in domain
         integer inumgarage     !number of garage buidligns in domain
         integer inumpolygon !MAN 08/10/2010 number of polygon nodes
! MAN 03/06/2007 fixed canopy variables type
         integer, allocatable :: canopy_ktop(:,:)
         real, allocatable :: canopy_top(:,:),canopy_atten(:,:,:)
         real, allocatable :: canopy_zo(:,:),canopy_ustar(:,:),canopy_d(:,:)
! end MAN 03/06/2007
         real umult_init,vmult_init !ANU vegetation canopy variables 08/04/05
         real ustar, lc, d, vk, ac !TMB 7/10/03 canopy variables
         real Az,Bz
         real uH,zd !ANU vegetation canopy variables 08/04/05

!         real, allocatable :: canopy_uin(:),canopy_uin_prof(:,:),umult_can(:), vmult_can(:)

!erp 1/3/2006
         real, allocatable :: atten(:)
!ERP 8/17/05 VARIABLES FOR building connectivity/ building group IDs
         integer, allocatable :: group_id(:)
!ERP end ID vars
!MAN 8/30/2005 stacked building fix
         real, allocatable :: zfo_actual(:)
         !integer, allocatable :: istart_canyon_N(:),iend_canyon_N(:),istart_canyon_S(:),iend_canyon_S(:) !MAN 7/5/2006
         !integer, allocatable :: jstart_canyon_E(:),jend_canyon_E(:),jstart_canyon_W(:),jend_canyon_W(:) !MAN 7/5/2006&
         !integer, allocatable :: kend_canyon_E(:),kend_canyon_N(:),kend_canyon_W(:),kend_canyon_S(:) !MAN 7/5/2006
! AAG 08/25/2006 Diffusion Variables
         integer diffusion_flag,diffstep,diffiter
         real dt                                       
		   real, allocatable :: visc(:,:,:), Fxd(:,:,:),Fyd(:,:,:) ,Fzd(:,:,:)
		   integer, allocatable :: rooftop_flag(:) ! AAG 09/13/06 added rooftop flag for individual building
		   real residual_reduction,domain_rotation,utmx,utmy
		   integer utmzone	
		   real, allocatable :: uo_roof(:,:,:),vo_roof(:,:,:) !NLB 09/12/05 For New Rooftop
! ERP parking garage additions
		   integer num_new_builds
! MAN negative buildings
         integer inumbuildneg
! MAN street canyon fix
         real max_velmag
! MAN 3/11/08 Staggered velocity write flag
         integer staggered_flag
! MAN 7/15/2008 Stretched vertical grid
         real, allocatable:: z(:),zm(:),dz_array(:)
         integer stretchgridflag
! MAN 8/17/2009 explosive source building damage
         integer damageflag
         integer, allocatable :: bld_damage(:)
         real explosionx,explosiony,hemass
         real Rdestroyed,Rdamaged,Rbuild
! MAN 8/25/2009 land use input
         integer landuse_flag,landuse_veg_flag,landuse_urb_flag,canopy_flag,lu_canopy_flag
         integer, allocatable :: landuse(:,:)
         real, allocatable :: landuse_height(:,:),landuse_atten(:,:)
! MAN 9/16/2010 polygon wake
         real, allocatable :: LrNode(:),LrFace(:)
! MAN 4/25/2013 polygon sidewall
         real, allocatable :: FaceRelWindDir(:)
         integer qcfd_flag
! MAN 2/10/2012 georeferenced sensor locations
         real meso_search_distance
         integer utmzoneletter
! MAN 3/06/2014 utc time
         integer utc_offset
         real dxy
      end module datamodule
