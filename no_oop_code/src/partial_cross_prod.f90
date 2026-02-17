subroutine pcrossprod (nrow,ncol,nrank, u, v, irow, pcol, nomega, r)

   implicit none

   ! arguments
   integer,          intent(in)  :: nrow, ncol, nrank, nomega
   integer,          intent(in)  :: irow(nomega), pcol(ncol+1)
   double precision, intent(in)  :: u(nrow,nrank), v(nrank,ncol)
   double precision, intent(out) :: r(nomega)

  ! variables
  integer :: j, ii, i


   do j = 1, ncol
      do ii = pcol(j)+1, pcol(j+1)          ! ii runs over the nnz of column j
         i        = irow(ii) + 1            ! row index is 0-based in R
         r(ii)    = dot_product(u(i,:), v(:,j))

      end do
   end do

end subroutine pcrossprod
