"""
Reads the column ids from a lammpstrj header
#### Returns:
- Int  column id for x
- Int  column id for y
- Int  column id for z
- Int  column id for element
"""
function get_headcolumns(head)
    if !("x" in keys(head["columns"]) )
        if "xu" in keys(head["columns"])
            head["columns"]["x"] = head["columns"]["xu"]
            head["columns"]["y"] = head["columns"]["yu"]
            head["columns"]["z"] = head["columns"]["zu"]
        end
    end
    xc, yc, zc = head["columns"]["x"], head["columns"]["y"], head["columns"]["z"]
    elc = head["columns"]["element"]
    return xc, yc, zc, elc
end

"""
Parse a lammpstrj file with elements and 3 coordinates per atom
#### Args:
- filename: String  path to lammpstrj file
##### Kwargs:
- start: Int  timestep to start reading from (default 0)
- stop: Int  timestep to stop reading at
- nth: Int  read every nth timestep (must be multiple of trajectory's timestep interval)
- ret_steps: Boolean  whether to return the timestepids as a array
- amsteps: Int  amount of steps to read (optional, will ignore stop)
#### Returns:
- elements: Array{String, 2}  elements[i, j] is the element of atom j in timestep i
- coords: Array{Float64, 3}  coords[i, j, k] is the kth coordinate of atom j in timestep i
- boxs: Array{Float64, 3}  boxs[i, j, k] is boxbounds in timestep i 
- (optional) steps: Array{Int, 1}  timestep ids

#### Notes:
The trajectory file has frames written at fixed intervals (delts). The nth parameter 
must be greater than or equal to this interval and must be a multiple of it.
"""
function parse_lammpstrj(
        filename::String; 
        start=0, stop=Inf, ret_steps=false, nth=nothing, amsteps=nothing)
    # read head of file (9 lines)
    head = Dict()
    open(filename) do file
        head = read_lammpstrj_head(file)
        firststep = head["timestep"]
        @debug "First Head contains: $(head)"
        if start == 0
            start = firststep
        end
        if head["timestep"] > start
            error("Start timestep is before first timestep in file")
        end

        delts, ts_byte = get_delts(file; ret_bytes=true)

        if nth == nothing
            nth = delts  # minimum amount of steps to skip, because inbetweens are not within the file
        else
            @assert nth >= delts "nth should be larger than the amount of steps inbeween the written steps to the file"
            @assert nth%delts == 0 "nth should be a multiple of the steps skipped written to file"
        end
        if stop != Inf
            @debug "Stopping at timestep $(stop)"
        else
            stop, last_pos = seekend_lammpstrj(file)
            ts_byte = last_pos÷((stop-firststep)÷delts-1)
            @debug "Set tsbyte by approx from final step to $(ts_byte)"
        end

                # skip over timesteps until start
        seekstart(file)
        if start > head["timestep"]
            start_byte = find_lammpstrj_timestep(file, start, approx_byte=ts_byte, delts=delts)
            @debug "found starting timestep at byte $(start_byte)"
        end
        @debug begin
            head = read_lammpstrj_head(file) 
            "Starting at timestep $(head["timestep"])"
        end
        @debug "file position: $(position(file))"

        if amsteps !== nothing
            @info "Using $amsteps instead of stop value"
            stop = start + (amsteps-1)*nth
        end
        timesteps = start:nth:stop
        @debug "reading the following timesteps: $(timesteps)"
        # initialize arrays to store data
        elements = Array{String}(undef, (stop-start)÷nth+1,  head["n_atoms"])
        coords = Array{Float64}(undef, (stop-start)÷nth+1,  head["n_atoms"], 3)
        boxs = Array{Float64}(undef, (stop-start)÷ nth+1,  2, 3)
        @debug "initialized arrays:", size(elements), size(coords), size(boxs)


        # read data
        @debug begin
            global starttime = time()
            "Reading total amount of $(length(timesteps)) timesteps"
        end
        if !("x" in keys(head["columns"]) )
            if "xu" in keys(head["columns"])
                @info "loading unwrapped coordinates"
            else
                error("No x column found in timestep")
            end
        end

        @debug "reading timesteps $(timesteps)"
        @debug "with the indices: $(1:(stop-start)÷(nth)+1)"
        for (tstep, tstepind) in zip(timesteps, 1:(stop-start)÷(nth)+1)
            if tstepind%1000 == 0
                @debug "Reading $(tstepind)th at time: $(time() - starttime)"
            end

            # TODO: swap with function call
            head = read_lammpstrj_head!(file)
            if head["timestep"] != tstep
                @debug "found tstep $(head["timestep"]) is not $(tstep)"
                @debug "Seeking to timestep $(tstep)"
                find_lammpstrj_timestep(file, tstep; approx_byte=ts_byte, delts=delts)
                head = read_lammpstrj_head!(file)
            end
            if !("x" in keys(head["columns"]) )
                xc, yc, zc = head["columns"]["xu"], head["columns"]["yu"], head["columns"]["zu"]
            else
                xc, yc, zc = head["columns"]["x"], head["columns"]["y"], head["columns"]["z"]
            end
            elc = head["columns"]["element"]
            for l in 1:head["n_atoms"]
                line = readline(file)
                line = split(line)  # split into array of Strings
                elements[tstepind, l] = line[elc]
                coords[tstepind, l, :] = [parse(Float64, x) for x in line[xc:zc]]
            end
            boxs[tstepind, :, :] = head["box_bounds"]

        end
        if ret_steps
            return elements, coords, boxs, [i for i in start:nth:stop] 
        end
        return elements, coords, boxs
    end
end

"""
Parse single timestep of lammpstrj file
#### Args:
- file: IOStream  file handle at any position
- tstep: Int  timestep to read
#### Kwargs:
- approx_byte: Int  approximate bytes per timestep 
- delts: Int  delta timestep per timestep written to lammpstrj file 
#### Returns:
- elements: Array{String, 1}  elements of atoms
- coords: Array{Float64, 2}  coordinates of atoms 
- box: Array{Float64, 2}  box bounds
"""
function parse_lammpstrj_step(file, tstep; approx_byte=nothing, delts=nothing)
    if eof(file)
        seekstart(file)
    end
    if isnothing(approx_byte) || isnothing(delts)
        delts, approx_byte = get_delts(file; ret_bytes=true)
    end
    find_lammpstrj_timestep(file, tstep; approx_byte=approx_byte, delts=delts)
    head = read_lammpstrj_head!(file)
    if ("x" in keys(head["columns"]) )
        xc, yc, zc = head["columns"]["x"], head["columns"]["y"], head["columns"]["z"]
    elseif ("xu" in keys(head["columns"]) )
        xc, yc, zc = head["columns"]["xu"], head["columns"]["yu"], head["columns"]["zu"]
    else
        @error "no x or xu columns found in file"
    end
    elc = head["columns"]["element"]
    elements = Array{String}(undef, head["n_atoms"])
    coords = Array{Float64}(undef, head["n_atoms"], 3)
    box = Array{Float64}(undef, 2, 3)

    for l in 1:head["n_atoms"]
        line = readline(file)
        line = split(line)  # split into array of Strings
        elements[l] = line[elc]
        coords[l, :] = [parse(Float64, x) for x in line[xc:zc]]
    end
    box[ :, :] = head["box_bounds"]
    return elements, coords, box
end

"""
Search for the last timestep and return the timestep number
#### Args:
- file: IOStream of a lammpstrj file at any position 
#### Returns:
- Int  timestep number of last timestep
- Int  position of last timestep in file
"""
function seekend_lammpstrj(file::IOStream)
    startpos = position(file)
    laststep, last_pos = seekend_lammpstrj!(file)
    seek(file, startpos)
    return laststep, last_pos   
end

function seekend_lammpstrj!(file::IOStream)
    seekend(file)
    @debug "Searching for final timestep"
    laststep = 0
    last_pos = 0
    s = ""
    while true
        skip(file, -2)
        c = read(file, Char)
        s = string(c, s)
        if c == '\n'
            s = ""
        end
        if endswith(s, "STEP")
            readline(file)
            last_pos = position(file) - 15  # position before "ITEM: TIMESTEP"
            laststep = parse(Int, readline(file))
            @debug "Final Timestep:", laststep
            break
        end
    end
    seek(file, last_pos)
    return laststep, last_pos   
end

"""
Read the head of a lammpstrj file timestep
#### Args:
- file: IOStream at position of head start
- offset: Int  offset from current position to start of head
Returns:
- Dict with keys:
  + timestep: Int
  + n_atoms: Int
  + box_bounds: Array{Float64, 2}  [xmin ymin zmin; xmax ymax zmax]
"""
function read_lammpstrj_head(file::IOStream; offset=0)
    pos = position(file)
    head = read_lammpstrj_head!(file; offset=offset)
    seek(file, pos)
    return head
end

"""
Read the head of a lammpstrj file timestep and update the file position
#### Args:
- file: IOStream at position of head start
##### Kwargs:
- offset: Int  offset from current position to start of head
#### Returns:
- Dict with keys:
  + timestep: Int
  + n_atoms: Int
  + box_bounds: Array{Float64, 2}  [xmin ymin zmin; xmax ymax zmax]
"""
function read_lammpstrj_head!(file::IOStream; offset=0)
    head = Dict()
    skip(file, offset)
    readline(file) # ITEM: TIMESTEP
    head["timestep"] = parse(Int, readline(file))
    readline(file) # ITEM: NUMBER OF ATOMS
    head["n_atoms"] = parse(Int, readline(file))
    readline(file) # ITEM: BOX BOUNDS
    head["box_bounds"] = [parse(Float64, x) for x in split(readline(file))]
    head["box_bounds"] = [head["box_bounds"] [parse(Float64, x) for x in split(readline(file))]]
    head["box_bounds"] = [head["box_bounds"] [parse(Float64, x) for x in split(readline(file))]]
    cols = split(readline(file))[3:end] # ITEM: ATOMS id type x y z
    head["columns"] = Dict()
    @debug "Columns in head: $(cols)"
    for (i, c) in enumerate(cols)
        head["columns"][c] = i
    end
    @debug "Read head: $(head)"
    return head
end

"""
Searches for start of defined Timestep in lammpstrj file
#### Args:
- file: IOStream at position of head start
- timestep: Int  timestep to search for
##### Kwargs:
- approx_byte: Int  approximate bytes per timestep 
- delts: Int  delta timestep per timestep written to lammpstrj file
- fsize: Int  size of file in bytes
- savety: Int  amount of timesteps to save in case of overshooting
Returns:
- position of timestep in file
"""
function find_lammpstrj_timestep(file::IOStream, timestep::Int; approx_byte=nothing, delts=nothing, fsize=nothing, savety=5)
    if eof(file)
        @info "File end found, skipping to start"
        seekstart(file)
    end
    pos = position(file)
    @debug "Starting search at byte position $(pos)"
    if approx_byte == nothing
        delts, approx_byte = get_delts(file; ret_bytes=true)
        @debug "Got delts and approxbyte"
        @debug delts, approx_byte
    elseif delts == nothing
        delts = get_delts(file)
        @debug delts
    end
    
    # Align to timestep boundary if needed
    while  !endswith(readline(file), "ITEM: TIMESTEP")
    end
    skip(file, -15)

    curr_head = read_lammpstrj_head(file)
    curr_ts = curr_head["timestep"]
    @debug "Current timestep: $(curr_ts)"
    skip_bytes = approx_byte*((timestep - curr_ts)÷delts-savety)
    seek_pos = position(file) + skip_bytes
    if seek_pos < 0
        seek_pos = 0
    end
    if fsize == nothing
        # get file size
        seekend(file)
        fsize = position(file)
    end
    seek(file, seek_pos)  # start approx 5 steps before position

    # check if we went to far
    if position(file) >= fsize-approx_byte
        laststep, pos = seekend_lammpstrj!(file)
        approx_byte = pos÷(laststep÷delts-1)
    end
    @debug "Searching for timestep $(timestep)"
    tofar = 0
    while !eof(file)
        l = readline(file)
        if endswith(l, "STEP")
            pos = position(file)
            step = parse(Int, readline(file))
            @debug "Found timestep $(step) at $(pos)"
            if step > timestep
                if tofar > 0
                    @error "Timestep $(timestep) is not in the file"
                    exit(1)
                end
                @debug "Approximation too far at $(step), seeking back"
                tofar = approx_byte*((step-timestep)÷delts+1)
                seek(file, pos - tofar)
                continue
            end
            if step == timestep
                seek(file, pos-length(l)-1)
                return pos - length(l) - 1
            end
        end
    end
    error("Did not find timestep in file")
end

"""
Determines the timestep interval (delts) between frames in a lammpstrj file
#### Args:
- file: IOStream  file handle at any position
- ret_bytes: Bool  whether to also return approximate bytes per timestep (default false)

#### Returns:
- Int  timestep interval between frames
- Int  (optional) approximate bytes per timestep if ret_bytes=true

#### Notes:
This function preserves the original file position.
"""
function get_delts(file::IOStream; ret_bytes::Bool=false)
    original_pos = position(file)
    
    # Find start of a timestep
    while !endswith(readline(file), "ITEM: TIMESTEP")
        if eof(file)
            seek(file, 0)  # Start from beginning if we hit EOF
        end
    end
    skip(file, -15)
    start_pos = position(file)
    
    # Read first timestep
    head1 = read_lammpstrj_head!(file)
    ts1 = head1["timestep"]
    
    # Skip atom lines of first timestep
    for _ in 1:head1["n_atoms"]
        readline(file)
    end
    
    # Calculate bytes for first timestep
    bytes_per_ts = position(file) - start_pos
    
    # Read second timestep if not at EOF
    if !eof(file)
        head2 = read_lammpstrj_head!(file)
        ts2 = head2["timestep"]
        delts = ts2 - ts1
    else
        # If we're at EOF, try from the beginning
        seek(file, 0)
        while !endswith(readline(file), "ITEM: TIMESTEP")
        end
        skip(file, -15)
        head1 = read_lammpstrj_head!(file)
        ts1 = head1["timestep"]
        
        # Skip atom lines
        for _ in 1:head1["n_atoms"]
            readline(file)
        end
        
        if !eof(file)
            head2 = read_lammpstrj_head!(file)
            ts2 = head2["timestep"]
            delts = ts2 - ts1
        else
            # If file has only one timestep, assume delts=1
            delts = 1
            @warn "Only one timestep found in file, assuming delts=1"
        end
    end
    
    # Restore original position
    seek(file, original_pos)
    
    return ret_bytes ? (delts, bytes_per_ts) : delts
end






