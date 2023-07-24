subroutine reorder_P(state_H, state_L, NumStates, RSteps, swaps) !input lowest state and it's higher swapped counterpart; the next state is the next state
integer NumStates, RSteps, i, j, swaps
double precision tmp(NumStates), P(NumStates, NumStates)

open(unit = 101, file = "P_Matrix.dat")
open(unit = 500, file = "reordered_P.dat")

do iR = RSteps, 1, -1
    read(101,*,END=100) tmp
    do i = 1, 80
        read(101,*,END=100) (P(i, j), j=1,NumStates)
    enddo

    tmp = P(:, state_L)
    P(:, state_L) = P(:, State_H)
    P(:, State_H) = tmp

    tmp = P(state_L, :)
    P(state_L, :) = P(State_H, :)
    P(state_H, :) = tmp

    tmp = P(:, state_L+1)
    P(:, state_L+1) = P(:, State_H-1)
    P(:, State_H-1) = tmp

    tmp = P(state_L+1, :)
    P(state_L+1, :) = P(State_H-1, :)
    P(state_H-1, :) = tmp

    do i = 1, NumStates
        do j = 1, NumStates
            tmp = P(i, State_L:(State_H-swaps))
            P(i,State_L:(State_L+1)) = P(i, (State_H-1):State_H)
            P(i, (State_L+swaps:NumStates)) = tmp
        enddo
    enddo
enddo
100 close(101)
close(500)

end subroutine reorder_P