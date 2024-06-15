program cluster_dist
implicit none

integer, parameter :: NumberOfParticles=1024, max_NumberOfParticles=1024, lx=200, ly=200,lz=200
integer :: i,j,k,nn,clustr,clustrr,clustr1,clustr2,no
integer :: clusteri
real*8, parameter :: llx=dfloat(lx), lly=dfloat(ly), llz=dfloat(lz)
integer, parameter :: n_iter = 5000, n_save = 50
integer, parameter :: n_total = n_iter/(2*n_save)
real*8, parameter :: dt_md = 0.01d0, time_step = 50.0d0
integer, parameter :: n_frame=40000, n_miss = 1!nint(time_step/dt_md)/n_save  
integer :: n_skip

!!real*8:: Up,U0,U1,U2,U3,U4

real*8 :: pos(3*NumberOfParticles),sigma_colloid=5.0d0!!,epsn,lQQ,lda,cut_hc2

!N_cluSter: number of clusters which contain i(how many) particles
!cluster_Index(i): non-zero means particle i belongs to at least one cluster
!cluster_partIndex(clustr,no_part(clustr)): store the particle Index for Cluster with Index "clustr",... 
!...and particle with Index "no_part(clustr)" in this cluster 
integer :: id_no, colorid(NumberOfParticles), norm_const(n_frame)
integer :: no_part(max_NumberOfParticles), cluster_partIndex(max_NumberOfParticles,max_NumberOfParticles)
integer :: cluster_Index(NumberOfParticles), N_cluSter(max_NumberOfParticles), bond_no
integer :: max_bond, bonds(NumberOfParticles), bond_dist(0:max_NumberOfParticles)
real*8 :: distx,disty,distz,dist2,lbond2
!!real*8 :: dp1x,dp1y,dp1z,dp2x,dp2y,dp2z,dp3x,dp3y,dp3z,dp4x,dp4y,dp4z,dp1,dp2,dp3,dp4
real*8 :: llxby2,llyby2,llzby2
integer :: clust_size(max_NumberOfParticles),total_clust_part
character :: element
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CHARACTER(*), PARAMETER :: ifileplace = "/Volumes/LaCie/ProteinDiffusion/DATA/Langevin_Q2D_SLAR_1024/A_phi_dot50/Epsilon2/"
CHARACTER(*), PARAMETER :: ofileplace = "/Volumes/LaCie/ProteinDiffusion/DATA/Langevin_Q2D_SLAR_1024/A_phi_dot50/Epsilon2/"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!lQQ=sigma_colloid*7.18d0

!!epsn=1.0d0
!!lda=1.794d0*sigma_colloid
!!cut_hc2=sigma_colloid*((2.0d0)**(0.002d0))
llxby2 = 0.5d0*llx; llyby2 = 0.5d0*lly; llzby2 = 0.5d0*llz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
lbond2 = (sigma_colloid*1.135795837d0)**2  !!!!some criterias for detecting cluster size
!(1.103708472,1.119680565,1.129086288,1.135795837)
!if (n_total/n_miss.lt.3000) then
!  n_frame = 5000!n_total/n_miss
! else if (n_total/n_miss.ge.3000) then
!  n_frame = 5000
!end if
! number of frames should be less than 3000
n_skip = 10000!n_total-n_frame*n_miss + (n_iter/n_save - n_total)

write(*,*) n_skip,n_total,n_frame,n_miss

  open (81, file=ofileplace//'cluster_dist.dat', status='unknown', form='formatted')
  open (82, file=ofileplace//'mov-com.xyz', status='unknown', form='formatted')
  open (83, file=ofileplace//'max_cluster_dist.dat', status='unknown', form='formatted')
  open (84, file=ofileplace//'bond_no.dat', status='unknown', form='formatted')! bonds/particle/frame
  open (85, file=ofileplace//'bond_dist.dat', status='unknown', form='formatted')! bond distribution at each frame
  open (86, file=ofileplace//'color_id.dat', status='unknown', form='formatted')
  open (87, file=ofileplace//'clust_size.dat', status='unknown', form='formatted')

!!Introducing a loop here to write many data?


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(079,file=ifileplace//'BDtraj101-StartStep-20000000.xyz')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) 'BDtraj32-data3.dat is under processing'
do j=1,n_skip
    read(79,*)    !skip the line "NumberOfParticles"
    read(79,*)    !skip the line "X Y Z" description
    read(79,*) pos!skip the first column and store the data to 'pos Array'
end do

clust_size=0
total_clust_part=0

do nn=1,n_frame
        do j=1,n_miss
  !  read(79,*)
    read(79,*) pos
        end do
  clustrr = 0
  clustr = 0
  no_part = 0
  cluster_Index = 0
  cluster_partIndex = 0
  bond_no = 0
  bonds = 0
  bond_dist = 0

  do i=1,NumberOfParticles-1
   do j=i+1,NumberOfParticles

      distx = pos(3*j-2) - pos(3*i-2)
      disty = pos(3*j-1) - pos(3*i-1)
      distz = pos(3*j) - pos(3*i)

      distx = distx - llx*anint(distx/llx)
      disty = disty - lly*anint(disty/lly)
      distz = distz - llz*anint(distz/llz)

      dist2 = distx*distx + disty*disty + distz*distz

    if (dist2.le.lbond2) then


        bond_no = bond_no + 2
        bonds(i) = bonds(i) + 1
        bonds(j) = bonds(j) + 1

        if (cluster_Index(i).eq.0.and.cluster_Index(j).eq.0) then

          clustrr = clustrr + 1  ! here we count it (only at the first time) as a cluster 
          clustr = clustrr   ! cluster (temporary) Index 

          no_part(clustr) = no_part(clustr) + 2 !accumulately count the Number of particles belongs to the same cluster

          cluster_partIndex(clustr,1) = i !
          cluster_partIndex(clustr,2) = j
          cluster_Index(i) = clustr !particle i or/and j belong to the same cluster with number cluster_Index
          cluster_Index(j) = clustr

        else if (cluster_Index(i).ne.0.and.cluster_Index(j).eq.0) then

          clustr = cluster_Index(i) ! passing the cluster index to the temporary variable clustr

          no_part(clustr) = no_part(clustr) + 1 !count the new particle j into the existing cluster "clustr", no_part +1

          cluster_partIndex(clustr,no_part(clustr)) = j !store the new particle index j into the existing cluster "clustr"
          cluster_Index(j) = clustr  !passing the cluster index to 

        else if (cluster_Index(i).eq.0.and.cluster_Index(j).ne.0) then

          clustr = cluster_Index(j)

          no_part(clustr) = no_part(clustr) + 1

          cluster_partIndex(clustr,no_part(clustr)) = i
          cluster_Index(i) = clustr

          !When two particles i and j already belong to two "different" clusters
          !(if they belong to the same cluster, nothing needs to be done)
          !Particle "i" belongs to cluster_Index(i) merges into cluster_Index(j) contains particle "j"
        else if (cluster_Index(i).ne.0.and.cluster_Index(j).ne.0) then
          !the particle index at the previous step
          clustr1 = cluster_Index(i)
          clustr2 = cluster_Index(j)

          if (clustr1.gt.clustr2) then

            clustr = clustr2  !the smaller cluster index  need

            do k=1,no_part(clustr1)  !loop all particles in "cluster_Index(i)" with bigger index

              no_part(clustr) = no_part(clustr) + 1
              !Take over all particles from "cluster_Index(i)", i.e., with bigger index
              cluster_partIndex(clustr,no_part(clustr)) = cluster_partIndex(clustr1,k) 
              cluster_Index(cluster_partIndex(clustr1,k)) = clustr!set these cluster index as "cluster_Index(j)"

            end do

            clustr = clustr1!what to do with the old cluster with bigger Index "cluster_Index(i)"?
            no_part(clustr) = 0! recount all the particles at cluster "cluster_Index(i)"

            do k=1,no_part(clustrr)

              no_part(clustr) = no_part(clustr) + 1 !recount all the particles at cluster "cluster_Index(i)"
              cluster_partIndex(clustr,no_part(clustr)) = cluster_partIndex(clustrr,k)
              cluster_Index(cluster_partIndex(clustrr,k)) = clustr

            end do

            no_part(clustrr) = 0  !the previous cluster is downsized into another cluster size
            clustrr = clustrr - 1

          else if (clustr2.gt.clustr1) then

            clustr = clustr1

            do k=1,no_part(clustr2)

              no_part(clustr) = no_part(clustr) + 1
              cluster_partIndex(clustr,no_part(clustr)) = cluster_partIndex(clustr2,k)
              cluster_Index(cluster_partIndex(clustr2,k)) = clustr

            end do

            clustr = clustr2
            no_part(clustr) = 0

            do k=1,no_part(clustrr)

              no_part(clustr) = no_part(clustr) + 1
              cluster_partIndex(clustr,no_part(clustr)) = cluster_partIndex(clustrr,k)
              cluster_Index(cluster_partIndex(clustrr,k)) = clustr

            end do

            no_part(clustrr) = 0
            clustrr = clustrr - 1

          end if
        end if
    endif
   end do
  end do
 !a frame has been checked
  N_cluSter = 0

  do clusteri=1,clustrr  !loop all clusters (at least 2 members), count number of clusters contain how many particles (no_part(i)=s) 

    N_cluSter(no_part(clusteri)) = N_cluSter(no_part(clusteri)) + 1

  end do

  do i=1,NumberOfParticles

    k = bonds(i) + 1!!!! represents number i has number of bonds(i). 
    !!!In distribution, the bonds number PLUS one, to leave a space to particles has zero bonds 
    bond_dist(k) = bond_dist(k) + 1!!!!

  end do

  max_bond = maxval(bonds) + 1

  do i=1,max_bond

    if (bond_dist(i).ne.0) then

      write(85,29) i-1,bond_dist(i)/dfloat(NumberOfParticles)
      29 format(2g18.10)

    end if
  end do

    write(82,24) NumberOfParticles
    24 format(1g18.10)

    write(82,25) 'description of xyz file'
    25 format(1g18.10)


  id_no = 0

   do k=1,NumberOfParticles

    if (cluster_Index(k).eq.0) then   !cluster_Index(k)=0 represents isolated particle doesn't form cluster

      write(82,21) '1',pos(3*k-2),pos(3*k-1),pos(3*k)
      21 format(4g18.10)

      id_no = id_no + 1   ! counting how many (also the index for colorid array) isolated particles
      colorid(id_no) = 1  ! colorid cluster size is 1 for this particle

    end if
   end do

    do i=1,clustrr !number of all different clusters
     do j=1,no_part(i)!for cluster 

      k = cluster_partIndex(i,j)!??????

      write(82,22) no_part(i),pos(3*k-2),pos(3*k-1),pos(3*k)  !particles with cluster size coded
      22 format(4g18.10)

      id_no = id_no + 1
      colorid(id_no) = no_part(i)   !colorid equals cluster size S 

     end do
    end do

    norm_const(nn) = maxval(colorid) ! largest cluster size of all frames

    write(86,31) (colorid-1.0d0)*0.7d0/9.0d0 !!!!!
    31 format(925g18.10)!!!925 numbers per row. Why???

    write(84,28) nn,dfloat(bond_no)/dfloat(NumberOfParticles) !averaging number of bonds per particle per frame, for percolation study?
    28 format(2g18.10)


   do i=1,max_NumberOfParticles

    if (N_cluSter(i).ne.0) then    !if N_cluSter(i) is larger than 0

      clust_size(i) = clust_size(i) + i*N_cluSter(i)
      total_clust_part = total_clust_part + i*N_cluSter(i) ! counting all the particles involved in clustering

      write(81,20) i,N_cluSter(i),i*N_cluSter(i)
      20 format(3g18.10)

    end if

    if (i.eq.maxval(no_part)) then

      write(83,30) nn*50.0d0,i     !!!!50 here means the frame jump. The maximum cluster size of each sampled frame
      30 format(3g18.10)

    end if
   end do

    write(81,26) !leave a line blank
    26 format(3g18.10)


end do!!!!end of reading the movie

write(*,*) 'come on'
  write(87,33) 1,(n_frame*NumberOfParticles-total_clust_part)/dfloat(NumberOfParticles*n_frame)
  33 format(2g18.10)

 do i=1,max_NumberOfParticles

  if (clust_size(i).ne.0) then

    write(87,32) i,clust_size(i)/dfloat(NumberOfParticles*n_frame)
    32 format(2g18.10)

  end if
 enddo

    write(*,*) maxval(norm_const)!!!!biggest cluster

  close(79)
  close(81)
  close(82)
  close(83)
  close(84)
  close(85)
  close(86)
  close(87)

end program cluster_dist
