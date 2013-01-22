      subroutine clique(sim,ds,left,rows,cols,alphac,top,bottom,orig)
c
      integer rows
      integer cols
      double precision sim(cols,cols)
      integer ds(rows,cols)
      integer left(rows)
      double precision alphac
      integer top
      integer bottom
      integer orig
c
      integer point
c    
c*************************** one ***************************
c
      top = 1
c
      do 10 i=1,cols-1
      where (sim(i,i+1:cols) < alphac) ds(top,i+1:cols) = 1
      if (all(ds(top,:) .eq. 1)) then
        cycle
      endif
      if (any(ds(top,:) .eq. 1)) then
        left(top) = i
        top = top + 1
      endif
   10 continue
c
      if (top .eq. 1) then
        return
      else
        orig = top - 1
      endif
c
c
c*************************** two **************************
c
      ds(top,left(1)) = 1
      ds(top,left(2)) = 1
c
      ds(top+1,:) = ds(2,:)
      ds(top+2,:) = ds(1,:)
      ds(top+3,:) = max(ds(1,:),ds(2,:))
      ds(top+1,left(1)) = 1
      ds(top+2,left(2)) = 1
c
      bottom = top + 3
      call collap(ds,top,bottom,rows,cols,orig)
c
      point = bottom 
c
      do 21 i=3,orig
        do 22 j=top,bottom
        if (point+2 > rows) then
          orig = -1
          return
        endif
        ds(point+1,:) = ds(j,:)
        ds(point+1,left(i)) = 1
        ds(point+2,:) = max(ds(i,:),ds(j,:))
        point = point + 2
   22   continue
      top = bottom + 1
      bottom = point
c     bottom = bottom * 2
      call collap(ds,top,bottom,rows,cols,orig)
      point = bottom
   21 continue
c
      end


      subroutine collap(ds,top,bottom,rows,cols,orig)
c
      integer rows
      integer cols
      integer ds(rows,cols)
      integer top 
      integer bottom 
      integer orig
c
      integer point
      integer keeper
c
c****************************** one *************************
c
      keeper = 0
      point = bottom 
c
      do 10 i=top,bottom
        do 11 j=top,i-1
        suma = sum(ds(i,:))
        sumb = sum(max(ds(i,:),ds(j,:)))
        if (all(ds(i,:)-ds(j,:) .eq. 0)) then
          goto 10
        else if (suma .eq. sumb) then
          goto 10
        endif
   11   continue
c
        do 13 j=i+1,bottom
        suma = sum(ds(i,:))
        sumb = sum(max(ds(i,:),ds(j,:))) 
        if (all(ds(i,:)-ds(j,:) .eq. 0)) then
          cycle
        else if (suma .eq. sumb) then
          goto 10
        endif
   13   continue
c
      point = point + 1
      ds(point,:) = ds(i,:)
      keeper = keeper + 1
   10 continue
c
c* collap ********************* two ******************************
c
      do 20 i=1,keeper
      ds(orig+i,:) = ds(bottom+i,:)
   20 continue
c
      top = orig + 1
      bottom = orig + keeper 
c
      return
c
      end
