subroutine inputpars
    use common_proc_icos
    coeff(1) ='ort0'
    coeff(2) ='newort'
    coeff(3) ='calc_randort_2f'
    coeff(4) ='calc_matrix_2f'
    coeff(5) ='calc_matrix_5f'
    coeff(6) ='transform_ort_2fTo5f'
    coeff(7) ='PR_threshold'
    coeff(8) ='calc_matrix_3f'
    coeff(9) ='transform_ort_2fTo3f'


    proc%calc_randort_2f='n'
    proc%calc_matrix_2f='n'
    proc%calc_matrix_5f='n'
    proc%transform_ort_2fTo5f='n'
    PR_threshold=180.

    return
end
