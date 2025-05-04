using HDF5, LinearAlgebra, Dates
using SINNLS


function load_nnls_array(filename::String)
    h5open(filename, "r") do file
        return read(file["A"]), read(file["B"]), read(file["dep_data"]), read(file["init_cols"]), read(file["w_pri"])
    end
end
         
function load_fail_array(filename::String)
    h5open(filename, "r") do file
        return read(file["A_fail"])
    end
end

function run_nnls(PREFIX::String, THREAD::Int)
    BLAS.set_num_threads(THREAD)

    A, B, dep_data, init_cols, w_pri = load_nnls_array("$PREFIX/matrix.h5")
    H = 3600.0


    n = size(A)[2]
    A_T = eltype(A)

    # Initialization
    x0_ = zeros(A_T, n)
    init_cols .+= 1 # 1-index for julia

    x0_[init_cols] = w_pri

    weights_final = nothing
    max_order = maximum(dep_data)

    for now_order in 0:max_order
        @info "$(Dates.now()) Now_order : $now_order"
        order_index = findall(x -> x <= now_order, dep_data)
        A_reduced = @view A[:, order_index]

        weights = vec(SI_NNLS(A_reduced, B,
                              x0_ = x0_[order_index],
                              total_time=min(16, 2.0 ^ (now_order - 1)) * H,
                              restart_ratio=0.7,
                              epi=1))

        if now_order < max_order
            x0_[order_index] = weights
        else
            weights_final = weights
        end
    end

    predict_suc_B = A * weights_final
    error = norm(predict_suc_B - B)
    b_norm = norm(B)

    A_fail = load_fail_array("$PREFIX/matrix.h5")
    predict_B = vcat(predict_suc_B, A_fail * weights_final)

    return error, b_norm, weights_final, predict_B
end