module Profiling

export init_profiling, track_memory, track_performance_start, track_performance_end, @time_operation, print_profiling_reports, cleanup_profiling, format_bytes, format_time

# Memory tracking variables - minimal overhead
const MEMORY_TRACKER = Ref{Dict{String, Int64}}(Dict{String, Int64}())
const PEAK_MEMORY = Ref{Int64}(0)
const INITIAL_MEMORY = Ref{Int64}(0)

# Performance tracking variables (basic stats)
const PERFORMANCE_TRACKER = Ref{Dict{String, Float64}}(Dict{String, Float64}())
const TOTAL_EXECUTION_TIME = Ref{Float64}(0.0)
const PERFORMANCE_START_TIME = Ref{Float64}(0.0)

# Simple profiling control flag
const PROFILING_ENABLED = Ref{Bool}(false)

function init_profiling(enabled::Bool=false)
    """
    Initialize simple profiling system.
    
    Arguments:
    - enabled: Enable basic timing and memory tracking
    """
    PROFILING_ENABLED[] = enabled
    
    if enabled
        PERFORMANCE_START_TIME[] = time()
        INITIAL_MEMORY[] = Base.gc_live_bytes()
        MEMORY_TRACKER[] = Dict{String, Int64}()
        PERFORMANCE_TRACKER[] = Dict{String, Float64}()
        PEAK_MEMORY[] = INITIAL_MEMORY[]
        
        println("Starting basic performance and memory tracking")
        println("Initial memory: $(format_bytes(INITIAL_MEMORY[]))")
    end
end

function track_memory(variable_name::String, size_bytes::Int64)
    """Track memory allocation for a specific variable with minimal overhead."""
    if !PROFILING_ENABLED[]
        return
    end
    
    MEMORY_TRACKER[][variable_name] = max(get(MEMORY_TRACKER[], variable_name, 0), size_bytes)
    current_memory = Base.gc_live_bytes()
    PEAK_MEMORY[] = max(PEAK_MEMORY[], current_memory)
end

function get_memory_size(var)
    """Get the memory size of a variable in bytes."""
    if !PROFILING_ENABLED[]
        return 0
    end
    
    try
        return Base.summarysize(var)
    catch
        return 0  # Fallback for variables that can't be sized
    end
end

function track_performance_start(operation_name::String)
    """Start timing a specific operation."""
    if !PROFILING_ENABLED[]
        return 0.0
    end
    return time()
end

function track_performance_end(operation_name::String, start_time::Float64)
    """End timing and record the duration for a specific operation."""
    if !PROFILING_ENABLED[]
        return 0.0
    end
    
    duration = time() - start_time
    PERFORMANCE_TRACKER[][operation_name] = get(PERFORMANCE_TRACKER[], operation_name, 0.0) + duration
    return duration
end

macro time_operation(operation_name, expr)
    """Simple macro for timing operations."""
    quote
        local result
        if PROFILING_ENABLED[]
            start_time = track_performance_start($(operation_name))
            result = $(esc(expr))
            track_performance_end($(operation_name), start_time)
        else
            result = $(esc(expr))
        end
        result
    end
end

function print_profiling_reports()
    """Print comprehensive profiling reports."""
    if !PROFILING_ENABLED[]
        return
    end
    
    # Calculate total execution time
    TOTAL_EXECUTION_TIME[] = time() - PERFORMANCE_START_TIME[]
    
    # Print performance and memory reports
    print_performance_report()
    print_memory_report()
    print_hamiltonian_construction_breakdown()
end

function print_performance_report()
    """Print a performance analysis report (top 10 operations)."""
    if !PROFILING_ENABLED[]
        return
    end
    
    println("\n" * "="^60)
    println("PERFORMANCE ANALYSIS REPORT")
    println("="^60)
    
    total_time = TOTAL_EXECUTION_TIME[]
    println("Total execution time: $(format_time(total_time))")
    
    if !isempty(PERFORMANCE_TRACKER[])
        println("\nTop 10 time-consuming operations:")
        println("-" ^ 50)
        
        # Filter out Hamiltonian sub-operations (keep only top-level operations)
        filtered_ops = Dict{String, Float64}()
        for (op_name, time_spent) in PERFORMANCE_TRACKER[]
            # Exclude detailed Hamiltonian sub-operations, keep main ones
            if !contains(lowercase(op_name), "hamiltonian_") || 
               op_name == "Hamiltonian_Construction_Total"
                filtered_ops[op_name] = time_spent
            end
        end
        
        # Sort filtered operations by time spent (descending) and take top 10
        sorted_ops = sort(collect(filtered_ops), by=x->x[2], rev=true)
        top_10_ops = sorted_ops[1:min(10, length(sorted_ops))]
        
        for (i, (operation, time_spent)) in enumerate(top_10_ops)
            percentage = (time_spent / total_time) * 100
            println("$(lpad(i, 2)). $(rpad(operation, 30)) $(format_time(time_spent)) ($(round(percentage, digits=1))%)")
        end
        
        # Show totals based on filtered operations (not double-counted)
        total_tracked_time = sum(values(filtered_ops))
        tracking_coverage = (total_tracked_time / total_time) * 100
        println("\nTotal tracked time: $(format_time(total_tracked_time)) ($(round(tracking_coverage, digits=1))% coverage)")
    end
    
    println("="^60)
end

function print_memory_report()
    """Print a memory usage report (top 10 consumers)."""
    if !PROFILING_ENABLED[]
        return
    end
    
    println("\n" * "="^60)
    println("MEMORY USAGE REPORT")
    println("="^60)
    
    # Overall memory statistics
    final_memory = Base.gc_live_bytes()
    memory_used = final_memory - INITIAL_MEMORY[]
    peak_memory = PEAK_MEMORY[]
    
    println("Initial memory: $(format_bytes(INITIAL_MEMORY[]))")
    println("Final memory: $(format_bytes(final_memory))")
    println("Memory used during execution: $(format_bytes(memory_used))")
    println("Peak memory usage: $(format_bytes(peak_memory))")
    
    # Top 10 memory consumers
    if !isempty(MEMORY_TRACKER[])
        println("\nTop 10 memory consumers:")
        println("-" ^ 40)
        
        # Filter out Hamiltonian sub-operations (keep only main operations)
        filtered_memory = Dict{String, Int64}()
        for (var_name, size_bytes) in MEMORY_TRACKER[]
            # Exclude detailed Hamiltonian sub-operations, keep main ones
            if !contains(lowercase(var_name), "hamiltonian_") || 
               contains(var_name, "Hamiltonian_Matrix") && !endswith(var_name, "_Total")
                filtered_memory[var_name] = size_bytes
            end
        end
        
        sorted_vars = sort(collect(filtered_memory), by=x->x[2], rev=true)
        top_10_vars = sorted_vars[1:min(10, length(sorted_vars))]
        
        for (i, (var_name, size_bytes)) in enumerate(top_10_vars)
            percentage = (size_bytes / peak_memory) * 100
            println("$(lpad(i, 2)). $(rpad(var_name, 25)) $(format_bytes(size_bytes)) ($(round(percentage, digits=1))%)")
        end
        
        # Show totals based on filtered memory (not double-counted)
        total_tracked = sum(values(filtered_memory))
        tracking_coverage = (total_tracked / peak_memory) * 100
        println("\nTotal tracked allocations: $(format_bytes(total_tracked)) ($(round(tracking_coverage, digits=1))% coverage)")
    end
    
    println("="^60)
end

function print_hamiltonian_construction_breakdown()
    """Print detailed breakdown of Hamiltonian construction steps with memory tracking."""
    if !PROFILING_ENABLED[]
        return
    end
    
    # Look for Hamiltonian-related operations
    hamiltonian_ops = Dict{String, Float64}()
    hamiltonian_matrix_memory = Dict{String, Int64}()
    hamiltonian_total_memory = Dict{String, Int64}()
    
    # Define the expected order of Hamiltonian construction steps
    hamiltonian_step_order = [
        "Hamiltonian_Matrix_Allocation",
        "Hamiltonian_Rotation_Terms", 
        "Hamiltonian_Electric_Field",
        "Hamiltonian_Magnetic_Field",
        "Hamiltonian_Spin_Rotation",
        "Hamiltonian_Electric_Dipolar",
        "Hamiltonian_Magnetic_Dipolar",
        "Hamiltonian_Hermitization"
    ]
    
    # Collect timing data
    for (op_name, time_spent) in PERFORMANCE_TRACKER[]
        if contains(lowercase(op_name), "hamiltonian")
            hamiltonian_ops[op_name] = time_spent
        end
    end
    
    # Collect memory data - separate matrix memory and total memory
    for (var_name, memory_size) in MEMORY_TRACKER[]
        if contains(lowercase(var_name), "hamiltonian")
            if endswith(var_name, "_Total")
                # Total memory usage - store with _Total suffix
                hamiltonian_total_memory[var_name] = memory_size
            else
                # Matrix memory only
                hamiltonian_matrix_memory[var_name] = memory_size
            end
        end
    end
    
    if !isempty(hamiltonian_ops)
        println("\n" * "="^120)
        println("HAMILTONIAN CONSTRUCTION BREAKDOWN")
        println("="^120)
        
        total_hamiltonian_time = 0.0
        hamiltonian_steps_found = String[]
        
        # Calculate total time for steps that actually ran
        for step_name in hamiltonian_step_order
            if haskey(hamiltonian_ops, step_name)
                total_hamiltonian_time += hamiltonian_ops[step_name]
                push!(hamiltonian_steps_found, step_name)
            end
        end
        
        total_execution_time = TOTAL_EXECUTION_TIME[]
        
        # Print header
        println("$(rpad("Step", 40)) $(rpad("Time", 12)) $(rpad("%", 8)) $(rpad("Matrix Size", 15)) $(rpad("Total Memory", 15))")
        println("-" ^ 120)
        
        # Print sub-steps with tree structure (only those that actually ran)
        step_counter = 0
        for step_name in hamiltonian_steps_found
            step_counter += 1
            time_spent = hamiltonian_ops[step_name]
            matrix_memory_used = get(hamiltonian_matrix_memory, step_name, 0)
            # Look for total memory with _Total suffix
            total_memory_used = get(hamiltonian_total_memory, step_name * "_Total", 0)
            
            percentage = total_hamiltonian_time > 0 ? (time_spent / total_hamiltonian_time) * 100 : 0.0
            
            # Format step name for display
            display_name = replace(step_name, "Hamiltonian_" => "")
            display_name = replace(display_name, "_" => " ")
            
            # Tree structure formatting
            tree_prefix = step_counter == length(hamiltonian_steps_found) ? "└── " : "├── "
            
            # Format output line
            time_str = format_time(time_spent)
            percent_str = "($(round(percentage, digits=1))%)"
            matrix_memory_str = matrix_memory_used > 0 ? format_bytes(matrix_memory_used) : "0 B"
            total_memory_str = total_memory_used > 0 ? format_bytes(total_memory_used) : "0 B"
            
            full_step_name = "     $(tree_prefix)$(display_name)"
            println("$(rpad(full_step_name, 40)) $(rpad(time_str, 12)) $(rpad(percent_str, 8)) $(rpad(matrix_memory_str, 15)) $(rpad(total_memory_str, 15))")
        end
        
        println("-" ^ 120)
        println("$(rpad("Total Hamiltonian construction:", 40)) $(rpad(format_time(total_hamiltonian_time), 12)) $(rpad("(100.0%)", 8))")
        
        # Memory summary
        if !isempty(hamiltonian_matrix_memory)
            println("\nHamiltonian memory summary:")
            final_matrix_memory = get(hamiltonian_matrix_memory, "Hamiltonian_Matrix_Final", 0)
            final_total_memory = get(hamiltonian_total_memory, "Hamiltonian_Matrix_Final_Total", 0)
            
            if final_matrix_memory > 0
                println("$(rpad("Final matrix size:", 40)) $(format_bytes(final_matrix_memory))")
            end
            if final_total_memory > 0
                println("$(rpad("Final total memory usage:", 40)) $(format_bytes(final_total_memory))")
            end
        end
        
        println("="^120)
    end
end

function cleanup_profiling()
    """Perform cleanup and final profiling tasks."""
    if !PROFILING_ENABLED[]
        return
    end
    
    # Simple cleanup with tracking
    if PEAK_MEMORY[] > 1024^3  # If peak memory > 1 GiB
        cleanup_start_time = time()
        println("\nSuggesting garbage collection for memory cleanup...")
        GC.gc()
        cleanup_time = time() - cleanup_start_time
        
        PERFORMANCE_TRACKER[]["Garbage_Collection"] = get(PERFORMANCE_TRACKER[], "Garbage_Collection", 0.0) + cleanup_time
        
        final_mem_after_gc = Base.gc_live_bytes()
        println("Memory after garbage collection: $(format_bytes(final_mem_after_gc))")
        println("Garbage collection took: $(format_time(cleanup_time))")
    end
end

function format_bytes(bytes::Int64)
    """Format bytes in human-readable format."""
    if bytes < 1024
        return "$(bytes) B"
    elseif bytes < 1024^2
        return "$(round(bytes/1024, digits=1)) KiB"
    elseif bytes < 1024^3
        return "$(round(bytes/1024^2, digits=1)) MiB"
    else
        return "$(round(bytes/1024^3, digits=2)) GiB"
    end
end

function format_time(seconds::Float64)
    """Format time in human-readable format."""
    if seconds < 1e-3
        return "$(round(seconds*1e6, digits=1)) μs"
    elseif seconds < 1.0
        return "$(round(seconds*1e3, digits=1)) ms"
    elseif seconds < 60.0
        return "$(round(seconds, digits=2)) s"
    elseif seconds < 3600.0
        minutes = floor(seconds / 60)
        secs = seconds - minutes * 60
        return "$(Int(minutes))m $(round(secs, digits=1))s"
    else
        hours = floor(seconds / 3600)
        minutes = floor((seconds - hours * 3600) / 60)
        return "$(Int(hours))h $(Int(minutes))m"
    end
end

end # module Profiling
