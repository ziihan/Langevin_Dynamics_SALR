program cluster_dist
implicit none

integer, parameter :: nofp=1560, max_npart=1560, lx=150, ly=150,lz=150
integer :: i,j,k,nn,clustr,clustrr,clustr1,clustr2,no
real*8, parameter :: llx=dfloat(lx), lly=dfloat(ly), llz=dfloat(lz)
integer, parameter :: n_iter = 2000000, n_save = 50
integer, parameter :: n_total = n_iter/(2*n_save)
real*8, parameter :: dt_md = 0.01d0, time_step = 50.0d0
integer, parameter :: n_miss = nint(time_step/dt_md)/n_save
integer :: n_frame, n_skip
real*8 :: pos(3*nofp),sigma_colloid=6.0d0
integer :: id_no, colorid(nofp), norm_const(500)
integer :: no_part(max_npart), cluster_part(max_npart,max_npart)
integer :: cluster_no(nofp), n_cluster(max_npart), bond_no
integer :: max_bond, bonds(nofp), bond_dist(0:max_npart)
real*8 :: distx,disty,distz,dist2,lbond2
real*8 :: llxby2,llyby2,llzby2
integer :: clust_size(max_npart),total_clust_part

llxby2 = 0.5d0*llx; llyby2 = 0.5d0*lly; llzby2 = 0.5d0*llz
lbond2 = (sigma_colloid*1.13268d0)**2

if (n_total/n_miss.lt.1000) then
  n_frame = n_total/n_miss
else if (n_total/n_miss.ge.1000) then
  n_frame = 1000
end if

n_skip = n_total-n_frame*n_miss + (n_iter/n_save - n_total)

write(*,*) n_skip,n_total,n_frame,n_miss


  open (81, file='cluster_dist.dat', status='unknown', form='formatted')
  open (82, file='mov.xyz', status='unknown', form='formatted')
  open (83, file='max_cluster_dist.dat', status='unknown', form='formatted')
  open (84, file='bond_no.dat', status='unknown', form='formatted')
  open (85, file='bond_dist.dat', status='unknown', form='formatted')
  open (86, file='color_id.dat', status='unknown', form='formatted')
  open (87, file='clust_size.dat', status='unknown', form='formatted')

  open(079,file='colloid_positions.dat')

do j=1,n_skip

    read(79,*)
    read(79,*) pos

end do

clust_size=0
total_clust_part=0

do nn=1,n_frame
  do j=1,n_miss-1

    read(79,*)
    read(79,*) pos

  end do

  clustrr = 0
  clustr = 0
  no_part = 0
  cluster_no = 0
  cluster_part = 0
  bond_no = 0
  bonds = 0
  bond_dist = 0

  do i=1,nofp-1
    do j=i+1,nofp

      distx = pos(3*j-2) - pos(3*i-2)
      disty = pos(3*j-1) - pos(3*i-1)
      distz = pos(3*j)   - pos(3*i)

      distx = distx - llx*anint(distx/llx)
      disty = disty - lly*anint(disty/lly)
      distz = distz - llz*anint(distz/llz)

      dist2 = distx*distx + disty*disty + distz*distz

      if (dist2.le.lbond2) then

        bond_no = bond_no + 2
        bonds(i) = bonds(i) + 1
        bonds(j) = bonds(j) + 1

        if (cluster_no(i).eq.0.and.cluster_no(j).eq.0) then

          clustrr = clustrr + 1
          clustr = clustrr
          no_part(clustr) = no_part(clustr) + 2
          cluster_part(clustr,1) = i
          cluster_part(clustr,2) = j
          cluster_no(i) = clustr
          cluster_no(j) = clustr

        else if (cluster_no(i).ne.0.and.cluster_no(j).eq.0) then

          clustr = cluster_no(i)
          no_part(clustr) = no_part(clustr) + 1
          cluster_part(clustr,no_part(clustr)) = j
          cluster_no(j) = clustr

        else if (cluster_no(i).eq.0.and.cluster_no(j).ne.0) then

          clustr = cluster_no(j)
          no_part(clustr) = no_part(clustr) + 1
          cluster_part(clustr,no_part(clustr)) = i
          cluster_no(i) = clustr

        else if (cluster_no(i).ne.0.and.cluster_no(j).ne.0) then

          clustr1 = cluster_no(i)
          clustr2 = cluster_no(j)

          if (clustr1.gt.clustr2) then

            clustr = clustr2

            do k=1,no_part(clustr1)

              no_part(clustr) = no_part(clustr) + 1
              cluster_part(clustr,no_part(clustr)) = cluster_part(clustr1,k)
              cluster_no(cluster_part(clustr1,k)) = clustr

            end do

            clustr = clustr1
            no_part(clustr) = 0

            do k=1,no_part(clustrr)

              no_part(clustr) = no_part(clustr) + 1
              cluster_part(clustr,no_part(clustr)) = cluster_part(clustrr,k)
              cluster_no(cluster_part(clustrr,k)) = clustr

            end do

            no_part(clustrr) = 0
            clustrr = clustrr - 1

          else if (clustr2.gt.clustr1) then

            clustr = clustr1

            do k=1,no_part(clustr2)

              no_part(clustr) = no_part(clustr) + 1
              cluster_part(clustr,no_part(clustr)) = cluster_part(clustr2,k)
              cluster_no(cluster_part(clustr2,k)) = clustr

            end do

            clustr = clustr2
            no_part(clustr) = 0

            do k=1,no_part(clustrr)

              no_part(clustr) = no_part(clustr) + 1
              cluster_part(clustr,no_part(clustr)) = cluster_part(clustrr,k)
              cluster_no(cluster_part(clustrr,k)) = clustr

            end do

            no_part(clustrr) = 0
            clustrr = clustrr - 1

          end if
        end if
      end if
    end do
  end do


  n_cluster = 0

  do i=1,clustrr

    n_cluster(no_part(i)) = n_cluster(no_part(i)) + 1

  end do

  do i=1,nofp

    k = bonds(i) + 1
    bond_dist(k) = bond_dist(k) + 1

  end do

  max_bond = maxval(bonds) + 1

  do i=1,max_bond

    if (bond_dist(i).ne.0) then

      write(85,29) i-1,bond_dist(i)/dfloat(nofp)
      29 format(2g18.10)

    end if
  end do

    write(82,24) nofp
    24 format(1g18.10)

    write(82,25) 'description of xyz file'
    25 format(1g18.10)


  id_no = 0

  do k=1,nofp

    if (cluster_no(k).eq.0) then

      write(82,21) '1',pos(3*k-2),pos(3*k-1),pos(3*k)
      21 format(4g18.10)

      id_no = id_no + 1
      colorid(id_no) = 1

    end if
  end do

  do i=1,clustrr
    do j=1,no_part(i)

      k = cluster_part(i,j)

      write(82,22) no_part(i),pos(3*k-2),pos(3*k-1),pos(3*k)
      22 format(4g18.10)

      id_no = id_no + 1
      colorid(id_no) = no_part(i)

    end do
  end do

    norm_const(nn) = maxval(colorid)

    write(86,31) (colorid-1.0d0)*0.7d0/9.0d0
    31 format(925g18.10)

    write(84,28) nn,dfloat(bond_no)/dfloat(nofp)
    28 format(2g18.10)


  do i=1,max_npart

    if (n_cluster(i).ne.0) then

      clust_size(i) = clust_size(i) + i*n_cluster(i)
      total_clust_part = total_clust_part + i*n_cluster(i)

      write(81,20) i,n_cluster(i),i*n_cluster(i)
      20 format(3g18.10)

    end if


    if (i.eq.maxval(no_part)) then

      write(83,30) nn*50.0d0,i
      30 format(3g18.10)

    end if
  end do

    write(81,26)
    26 format(3g18.10)

end do

  write(87,33) 1,(n_frame*nofp-total_clust_part)/dfloat(nofp*n_frame)
  33 format(2g18.10)

do i=1,max_npart

  if (clust_size(i).ne.0) then

    write(87,32) i,clust_size(i)/dfloat(nofp*n_frame)
    32 format(2g18.10)

  end if
enddo

write(*,*) maxval(norm_const)

  close(79)
  close(81)
  close(82)
  close(83)
  close(84)
  close(85)
  close(86)
  close(87)

end program cluster_dist

