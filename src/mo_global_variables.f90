module mo_global_variables
  !
  use mo_kind, only: i4, dp

  implicit none

  public

  ! date definition for starting and end dates
  type period
    integer(i4) :: y_start, m_start, d_start
    integer(i4) :: y_end, m_end, d_end
    integer(i4) :: j_start, j_end
    integer(i4) :: n_years, n_months, n_days
    integer(i4) , allocatable :: time_points(:)
  end type period

contains

  ! initialize object date and calculate derived variables from:
  ! y_start, m_start, d_start, y_end, m_end, d_end
  type (period) function period_init(y_start, m_start, d_start, y_end, m_end, d_end, time_points)

    use mo_julian, only: date2dec
    
    implicit none
    
    integer(i4), intent(in) :: y_start, m_start, d_start
    integer(i4), intent(in) :: y_end, m_end, d_end
    integer(i4), intent(in) :: time_points(:)

    period_init%y_start = y_start
    period_init%m_start = m_start
    period_init%d_start = d_start
    period_init%y_end = y_end
    period_init%m_end = m_end
    period_init%d_end = d_end
    period_init%time_points = time_points

    period_init%j_start = date2dec(period_init%d_start, period_init%m_start, period_init%y_start)
    period_init%j_end = date2dec(period_init%d_end, period_init%m_end, period_init%y_end)

    period_init%n_years = period_init%y_end - period_init%y_start + 1_i4
    period_init%n_days = period_init%j_end - period_init%j_start + 1_i4
    period_init%n_months = (12_i4 - period_init%m_start + 1_i4) &
        + (period_init%y_end - period_init%y_start - 1_i4) * 12_i4 &
        + period_init%m_end
    
  end function period_init
  
end module mo_global_variables
