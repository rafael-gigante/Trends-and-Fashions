program trends_and_fashions
    integer(8), parameter :: N = 10 ** 5, L = 10 ** 3, pcrit =  10 ** 1
    integer, dimension(N) :: Li 
    integer, dimension(L) :: N1 = 0, N2 = 0, p
    real(8) :: S, Pc, Sc
    real(8), dimension(L) :: ni
    real(8), dimension(800) :: function_data, period

    call srand(53)

    !Abertura de arquivos
    open(1, file='data.out')
    open(2, file='entropy.out')
    open(3, file='percolation.out')
    open(4, file='susceptibility.out')
    open(5, file='label_decay.out')
    open(6, file='period.out')

    !Gerando a condição inicial do sistema
    do i = 1, N
        Li(i) = int(L * rand()) + 1
        N1(Li(i)) = N1(Li(i)) + 1
        Li(i) = int(L * rand()) + 1
        N2(Li(i)) = N2(Li(i)) + 1
    end do
    do i = 1, L
        p(i) = N2(i) - N1(i)
    end do
    write(1,*)
    write(1,*)
    N1 = N2

    icount = 0
    !Dinâmica
    do it = 1, 800
        do k = 1, N
            i = int(N * rand() + 1)
            j = int(N * rand() + 1)
            if (p(Li(i)) < p(Li(j))) then
                N2(Li(i)) = N2(Li(i)) - 1
                N2(Li(j)) = N2(Li(j)) + 1
                Li(i) = Li(j)
            else if (p(Li(i)) < pcrit) then
                m = int(L * rand() + 1)
                if (N2(m) == 0) then
                    N2(Li(i)) = N2(Li(i)) - 1             
                    N2(m) = N2(m) + 1
                    Li(i) = m
                end if
            end if 
        end do

        Pc = maxval(N2)
        function_data(it) = Pc / N
        S = 0
        Sc = 0
        do i = 1, L
            p(i) = N2(i) - N1(i)
            ni(i) = (real(N2(i)) / real(N))
            S = S -  ni(i) * log(ni(i))
            Sc = Sc + (N2(i) ** 2)
            if (it >= 200) then
                write(1,*) i, N2(i)
            end if
        end do
        Sc = Sc - Pc ** 2
        
        if (maxval(ni) == 1 .and. icount == 0) then
            icount = it
        end if 
        write(1,*)
        write(1,*)
        write(2,*) it, S
        write(3,*) it, (Pc / N)
        write(4,*) it, Sc / N ** 2
        if (icount /= 0) then
            write(5,*) (it - icount), maxval(ni)
        end if
        N1 = N2
    end do

    peak_count = 0
    sum_periods = 0.0
    
    ! Encontrar os picos na função
    do i = 2, 800 - 1
        if (function_data(i) > function_data(i - 1) .and. function_data(i) > function_data(i + 1)) then
            peak_count = peak_count + 1
            period(int(peak_count)) = i
        else if ((function_data(i) == function_data(i-1)) .and. function_data(i) > function_data(i+1)) then
            peak_count = peak_count + 1
            period(int(peak_count)) = i
        end if
        
    end do

    ! Calcular a média dos períodos entre os picos consecutivos
    if (peak_count > 1) then
        do i = 1, int(peak_count) - 1
            sum_periods = sum_periods + (period(i + 1) - period(i))
            write(6,*) i, (period(i + 1) - period(i))
        end do
        peak_period = sum_periods / real(peak_count - 1)
        print *, "O período aproximado da função é:", peak_period
    else
        print *, "Não foi possível encontrar picos suficientes para calcular o período."
    end if


end program trends_and_fashions