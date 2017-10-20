"""
    pushmetadata(opt, key, val)

Push the given metadata chunk `opt`, which is a dictionary of vectors of
metadata chunks.
"""
function pushmetadata(opt, key, val)
    symbkey = Symbol(strip(String(key)))
    if symbkey in keys(opt)
        push!(opt[symbkey], val)
    else
        opt[symbkey] = [val]
    end
end