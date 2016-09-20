      subroutine permute(class,pclass,numplt,tclass)
c
      integer numplt
      integer class(numplt)
      integer pclass(numplt)
c
c* local
c
      integer tclass(numplt)
      integer pool
      integer index
      double precision unifrnd
c
      call rndstart()
c
      pool = numplt
      do 10 i=1,numplt
      tclass(i) = class(i)
   10 continue
c
      do 11 i=1,numplt
      index = ceiling(unifrnd()*(pool))
      pclass(i) = tclass(index)
      tclass(index) = tclass(pool)
      pool = pool - 1
   11 continue
c
      call rndend
c
      return
c
      end

