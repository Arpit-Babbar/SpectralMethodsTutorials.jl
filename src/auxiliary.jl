using JSON3

src_dir() = @__DIR__ # Directory of file
test_dir() = src_dir() * "/../test"
test_data_dir() = test_dir() * "/data"

# convert struct to a named tuple
to_named_tuple(p) = (; (v=>getfield(p, v) for v in fieldnames(typeof(p)))...)

# convert named tuple to dict
to_dict(nt::NamedTuple) = Dict(string.(keys(nt)) .=> values(nt))

# convert struct to dict
to_dict(p::Any) = to_dict(to_named_tuple(p))

# write dictionary as JSON
function write_dict(d::Dict, filename::String)
   dict_file = open(filename, "w")
   JSON3.write(dict_file, d)
   close(dict_file)
   return nothing
end

# read dictionary from JSON format
function read_dict(filename::String)
   dict_file = open(filename)
   dict_data = JSON3.read(dict_file)
   close(dict_file)
   return dict_data
end

# write struct
function write_struct(obj::Any, filename::String)
   dict = to_dict(obj)
   write_dict(dict, filename)
   return nothing
end

# compare dictionaries
function compare_dict(d1_::Union{Dict, JSON3.Array}, d2_::Any)
   d1 = JSON3.read(JSON3.write(d1_))
   d2 = JSON3.read(JSON3.write(d2_))
   test_result = true
   for key in keys(d1)
      if !(d1[key] ≈ d2[key])
         test_result = false
      end
   end
   return test_result
end

function test_object_with_data(test, obj::Any, data_name::String)
   data_loc = test_data_dir() * "/" * data_name
   if test == false
      write_struct(obj, data_loc)
   end
   data = to_dict(obj)
   expected_data = read_dict(data_loc)
   test_result = compare_dict(data, expected_data)
   return test_result
end

# TODO - Write a unit test just for this function
function flatten_indices(i,j,N)
   n = i + (j-1)*(N-1)
   return n
end
