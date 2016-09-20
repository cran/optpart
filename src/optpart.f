      subroutine optpart(simptp,simptc,simctc,simrat,cardin,
     +                  numplt,numclu,musubx,final,clusid,
     +                  maxdmu,maxitr,numitr)
c
c* passed in/out
c
      integer numplt                            ! number of plots
      integer numclu                            ! number of clusters
      integer maxitr
      integer numitr
      double precision simptp(numplt,numplt)    ! similarity plot-to-plot
      double precision simptc(numplt,numclu)    ! similarity of plot-to-cluster
      double precision simctc(numclu,numclu)    ! similarity of cluster-to-cluster
      double precision simrat(maxitr)           ! overal similarity ratio
      double precision cardin(numclu)           ! total membership of cluster
      double precision musubx(numplt,numclu)    ! membership of plot in cluster
      double precision final(numplt,numclu)
      integer clusid(numplt)                    ! cluster ID for plot
      double precision maxdmu                   ! maximum membership transfer
c
c* common numer 
c
      double precision sumswi
      double precision sumsbt
      double precision sumcwi
      double precision sumcbt
c
      common / numer /  sumswi,sumsbt,sumcwi,sumcbt
c  
      character rorc*1
c
c* local
c
      double precision dratio
      double precision maxsim
      double precision tmprat
c
c**************************** one *********************************
c
      do 10 i=1,numplt
      if (clusid(i) .gt. 0) then
        musubx(i,clusid(i)) = 1.0
        do 11 j=1,numclu
        if (j .eq. clusid(i)) goto 11
        musubx(i,j) = 0.0
   11   continue
      endif
   10 continue
c
      rorc = 'R'
c
c* npmmds ************************* three ******************************
c
      do 40 i=1,maxitr         
   45 call calcar(musubx,cardin,numplt,numclu)
      call fclctc(simptp,simctc,numplt,numclu,musubx)
      call ratio(simctc,numclu,cardin,tmprat)
      simrat(i) = tmprat
      if (i .eq. 1) then
        dratio = 0.0
      else
        dratio = simrat(i) - simrat(i-1)
      endif
c
      if (dratio .gt. 0.0 .or. i .eq. 1) then
        do 43 j=1,numplt
          do 48 k=1,numclu
          final(j,k) = musubx(j,k)
   48     continue
   43   continue
      else if (dratio .le. 0 .and. i .gt. 1) then
        if (rorc .eq. 'R') then
          rorc = 'C'
        else if (rorc .eq. 'C') then
          rorc = 'U'
        else if (rorc .eq. 'U') then
          numitr = i
          goto 41
        endif
        do 46 j=1,numplt
          do 47 k=1,numclu
          musubx(j,k) = final(j,k)
   47     continue
   46   continue
        goto 45
      endif
c
      numitr = maxitr
c
      call fclptc(simptp,simptc,numplt,numclu,musubx,cardin)
      call deltam(simptc,musubx,numplt,numclu,maxdmu,cardin,rorc)
   40 continue
c
   41 continue
c
c* altfcl ****************************** five *******************************
c
      call fclctc(simptp,simctc,numplt,numclu,musubx)
c
      do 50 i=1,numplt
      maxsim = 0.0
        do 51 j=1,numclu
        if (musubx(i,j) .gt. maxsim) clusid(i) = j
   51   continue
   50 continue
c
      end
c
c* dsv_fcl ************** subroutine fclctc *************************
c
      subroutine fclctc(simptp,simctc,numplt,numclu,musubx)
c
c* passed
c
      integer numplt
      integer numclu
      double precision simptp(numplt,numplt)
      double precision simctc(numclu,numclu)
      double precision musubx(numplt,numclu)
c
c* local
c
      double precision weight
      double precision sumnum
      double precision sumden
c
c* dsv_fcl/fclctc *************** one **********************************
c
      do 10 i=1,numclu
        do 11 j=i,numclu
        simctc(i,j) = 0.0
        sumnum = 0.0
        sumden = 0.0
          do 12 k=1,numplt
          if (musubx(k,i) .le. 0.0) goto 12 
            do 13 l=1,numplt
            if (k .eq. l) goto 13
            if (musubx(l,j) .le. 0.0) goto 13
            weight = min(musubx(k,i),musubx(l,j))
            sumnum = sumnum + (simptp(k,l) * weight)
            sumden = sumden + weight
   13       continue
   12     continue
        if (sumden .lt. 0.01) then
          simctc(i,j) = 0.0
        else
          simctc(i,j) = sumnum /sumden
        endif
        simctc(j,i) = simctc(i,j)
   11   continue
   10 continue
c
      return
c
      end
c
c* dsv_fcl ******************** subroutine fclptc ***************************
c
      subroutine fclptc(simptp,simptc,numplt,numclu,musubx,cardin)
c
c* passed
c
      integer numplt
      integer numclu
      double precision simptp(numplt,numplt)
      double precision simptc(numplt,numclu)
      double precision musubx(numplt,numclu)
      double precision cardin(numclu)
c
c* local
c
      double precision sumnum
      double precision sumden
c
c* dsv_fcl/fclptc ****************** one ************************************
c
      do 10 i=1,numplt
        do 11 j=1,numclu
        if (cardin(j) .eq. 0.0) then
          simptc(i,j) = 0.0
        else
          sumnum = 0.0
          sumden = 0.0
            do 12 k=1,numplt
            if (musubx(k,j) .le. 0.0) goto 12
            if (i .eq. k) goto 12
            sumnum = sumnum + simptp(i,k) * musubx(k,j)
            sumden = sumden + musubx(k,j)
   12       continue
          if (sumden .le. 0.0) sumden = 1.0
          simptc(i,j) = sumnum / sumden
        endif
   11   continue
   10 continue
c
      return

      end
c
c* dsv_fcl ********* subtoutine ratio *********************************
c
      subroutine ratio(simctc,numclu,cardin,simrat)
c
c* passed
c
      integer numclu
      double precision simctc(numclu,numclu)
      double precision cardin(numclu)
      double precision simrat
c
c* common numer
c
      double precision sumswi
      double precision sumsbt
      double precision sumcwi
      double precision sumcbt
c
      common / numer /  sumswi,sumsbt,sumcwi,sumcbt
c
c
c* local
c
      double precision weight
      double precision simwis
      double precision simbts
c
c* dsv_fcl/ratio ************** one **********************************
c          
      sumswi = 0.0
      sumsbt = 0.0
      sumcwi = 0.0
      sumcbt = 0.0
c
      do 10 i=1,numclu
      if (cardin(i) .eq. 0.0) goto 10
        do 11 j=i,numclu
        if (cardin(j) .eq. 0.0) goto 11
        if (i .eq. j) then
          weight = (cardin(i)**2-cardin(i))/2
          sumswi = sumswi + weight * simctc(i,i)
          sumcwi = sumcwi + weight
        else
          weight = cardin(i) * cardin(j)
          sumsbt = sumsbt + weight * simctc(i,j)
          sumcbt = sumcbt + weight
        endif
   11   continue
   10 continue
c
      simwis = sumswi / sumcwi
      simbts = sumsbt / sumcbt
      simrat = simwis / simbts
c
      return
c
      end
c
c* dsv_fcl ************ subroutine deltam ****************************
c
      subroutine deltam(simptc,musubx,numplt,numclu,maxdmu,cardin,rorc)
c
      parameter (maxclu=100)
c
c* passed
c
      integer numplt
      integer numclu
      double precision simptc(numplt,numclu)
      double precision musubx(numplt,numclu)
      double precision maxdmu
      double precision cardin(numclu)
      character rorc*1
c
c* common numer
c
      double precision sumswi
      double precision sumsbt
      double precision sumcwi
      double precision sumcbt
c
      common / numer /  sumswi,sumsbt,sumcwi,sumcbt
c
c* local
c
      double precision dmusub
      double precision newswi
      double precision newsbt
      double precision newcwi
      double precision newcbt
      double precision ratio
      double precision musj,musk
      double precision dratio(maxclu,maxclu)
c
      double precision maxdif(maxclu,maxclu)
      double precision dltamu(maxclu,maxclu)
      double precision tmpmax
      integer pltdif(maxclu,maxclu)
      integer pltnum
      double precision summu
c
c* dsv_fcl/deltam ********** one **********************************
c
      ratio = (sumswi/sumcwi) / (sumsbt/sumcbt)
c
      do 18 i=1,numclu
        do 19 j=1,numclu
        maxdif(i,j) = 0.0
        dltamu(i,j) = 0.0
        pltdif(i,j) = 0
   19   continue
   18 continue
c
      do 10 i=1,numplt
        do 11 j=1,numclu
        musj = musubx(i,j)
        if (musj .le. 0) goto 11
          do 12 k=1,numclu
          dratio(j,k) = 0.0
          musk = musubx(i,k)
          if (j .eq. k) goto 12
          if (musk .lt. 0) goto 12
          dmusub = min(musj,1.0-musk) * maxdmu
          if (dmusub .gt. 0) then
            newswi = sumswi-((simptc(i,j)*(cardin(j)-musj))*dmusub)
     +                     +((simptc(i,k)*(cardin(k)-musk))*dmusub)
            newcwi = sumcwi - ((cardin(j)-musj)*dmusub)
     +                      + ((cardin(k)-musk)*dmusub)
            newsbt = sumsbt+((simptc(i,j)*(cardin(j)-musj))*dmusub)
     +                     -((simptc(i,k)*(cardin(k)-musk))*dmusub)
            newcbt = sumcbt + ((cardin(j)-musj)*dmusub)
     +                      - ((cardin(k)-musk)*dmusub)
            dratio(j,k) = ((newswi/newcwi) / (newsbt/newcbt)) /
     +                ratio
            if (dratio(j,k) .gt. maxdif(j,k)) then
              maxdif(j,k) = dratio(j,k)
              pltdif(j,k) = i
              dltamu(j,k) = dmusub
            endif
          endif
   12     continue
   11   continue
   10 continue
c
c* dsv_fcl/deltam ************** two **********************************
c
      tmpmax = 1.0
      numrow = 0
      numcol = 0
c
   20 do 21 i=1,numclu
        do 22 j=1,numclu
        if (i .eq. j) goto 22
          if (maxdif(i,j) .gt. tmpmax) then
            tmpmax = maxdif(i,j)
            numrow = i
            numcol = j
          endif
   22     continue
   21   continue
c
        if (tmpmax .gt. 1.0) then
          pltnum = pltdif(numrow,numcol)
          dmusub = min(musubx(pltnum,numrow),1.0-musubx(pltnum,numcol),
     +        dltamu(numrow,numcol))
     +           * maxdmu
          musubx(pltnum,numrow) =
     +    musubx(pltnum,numrow) - dmusub
          musubx(pltnum,numcol) =
     +        musubx(pltnum,numcol) + dmusub
          if (rorc .eq. 'C') then
            do 23 i=1,numclu
            maxdif(numrow,i) = 0.0
            maxdif(i,numrow) = 0.0
            maxdif(i,numcol) = 0.0
            maxdif(numcol,i) = 0.0
   23      continue
          else if (rorc .eq. 'U') then
            tmpmax = 1.0
            goto 26
          else
            maxdif(numrow,numcol) = 0.0
            maxdif(numcol,numrow) = 0.0
          endif
          tmpmax = 1.0
          goto 20
        endif
c
   26   do 24 i=1,numplt
        summu = 0.0
          do 25 j=1,numclu
          summu = summu + musubx(i,j)
   25     continue
   24   continue
c
        return
c
        end
c
c* dsv_fcl ******************** subroutine calcar ******************
c
      subroutine calcar(musubx,cardin,numplt,numclu)
c
c* passed
c
      integer numplt
      integer numclu
      double precision musubx(numplt,numclu)
      double precision cardin(numclu)
c
c* local
c
      double precision total
c
c* dsv_fcl/calcar ************** one ******************************
c
      do 10 i=1,numclu
      cardin(i) = 0.0
   10 continue
c
      do 11 i=1,numplt
        do 12 j=1,numclu
        cardin(j) = cardin(j) + musubx(i,j)
        total = total + musubx(i,j)
   12   continue
   11 continue
c
      return
c
      end

