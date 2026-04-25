import os

def list_filenames_in_folders():
    """
    Recursively list all filenames in specified folders and save to a file
    """
    
    # Define the step folders
    step_input_folder_name = []
    step_input_folder_name.append("Step1")
    step_input_folder_name.append("Step2")
    step_input_folder_name.append("Step3")
    step_input_folder_name.append("Step4_NEW")
    step_input_folder_name.append("Step5_New")
    step_input_folder_name.append("Step6_New")
    
    # Root folder where script resides (current directory)
    root_folder = os.getcwd()
    
    # Output file path
    output_file = os.path.join(root_folder, "all_filenames.txt")
    
    # Also create a detailed CSV file with folder information
    detailed_output_file = os.path.join(root_folder, "all_filenames_detailed.csv")
    
    print(f"Script running from: {root_folder}")
    print(f"Output will be saved to: {output_file}")
    print(f"Detailed output will be saved to: {detailed_output_file}")
    print("-" * 60)
    
    # Store all filenames
    all_files = []
    detailed_files = []
    
    # Track statistics
    total_files = 0
    folders_processed = 0
    
    # Process each step folder
    for step_idx, step_name in enumerate(step_input_folder_name, start=1):
        # Construct the full path
        folder_path = os.path.join("data", "cluster", "endsim", step_name)
        absolute_path = os.path.join(root_folder, folder_path)
        
        print(f"Step {step_idx}: {step_name}")
        print(f"  Looking in: {absolute_path}")
        
        # Check if folder exists
        if not os.path.exists(absolute_path):
            print(f"  WARNING: Folder does not exist! Skipping...")
            print(f"  Expected path: {absolute_path}")
            continue
        
        folders_processed += 1
        files_in_this_folder = []
        
        # Walk through all subdirectories and files
        for dirpath, dirnames, filenames in os.walk(absolute_path):
            # Get relative path from the step folder
            rel_path = os.path.relpath(dirpath, absolute_path)
            if rel_path == '.':
                rel_path = step_name
            
            for filename in filenames:
                full_path = os.path.join(dirpath, filename)
                file_size = os.path.getsize(full_path)
                total_files += 1
                files_in_this_folder.append(filename)
                
                # For simple list (just filenames)
                all_files.append(filename)
                
                # For detailed list (with path and metadata)
                detailed_files.append({
                    'step': step_name,
                    'step_number': step_idx,
                    'file_path': os.path.join(folder_path, rel_path, filename) if rel_path != step_name else os.path.join(folder_path, filename),
                    'filename': filename,
                    'full_absolute_path': full_path,
                    'size_bytes': file_size,
                    'subdirectory': rel_path if rel_path != step_name else ''
                })
        
        print(f"  Found {len(files_in_this_folder)} files")
        print(f"  Total files so far: {total_files}")
        print()
    
    # Write simple list of filenames (with duplicates if same filename appears in multiple folders)
    with open(output_file, 'w') as f:
        f.write("# List of all filenames found in Step1-Step6 folders\n")
        f.write("# Note: Same filename may appear multiple times if duplicated across folders\n")
        f.write("# For detailed information with paths, see all_filenames_detailed.csv\n")
        f.write("#" * 60 + "\n\n")
        
        for filename in sorted(set(all_files)):  # Use set() to get unique filenames
            f.write(f"{filename}\n")
    
    # Write detailed CSV file with full information
    import csv
    with open(detailed_output_file, 'w', newline='') as csvfile:
        if detailed_files:
            fieldnames = ['step_number', 'step', 'subdirectory', 'filename', 'file_path', 'full_absolute_path', 'size_bytes']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            # Sort by step number then filename
            detailed_files_sorted = sorted(detailed_files, key=lambda x: (x['step_number'], x['filename']))
            for file_info in detailed_files_sorted:
                writer.writerow(file_info)
    
    # Also create a summary text file
    summary_file = os.path.join(root_folder, "file_summary.txt")
    with open(summary_file, 'w') as f:
        f.write("FILE SCAN SUMMARY\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Scan performed from: {root_folder}\n")
        f.write(f"Total folders processed: {folders_processed}/6\n")
        f.write(f"Total files found: {total_files}\n")
        f.write(f"Unique filenames: {len(set(all_files))}\n\n")
        
        f.write("Files per folder:\n")
        for step_idx, step_name in enumerate(step_input_folder_name, start=1):
            folder_path = os.path.join("data", "cluster", "endsim", step_name)
            absolute_path = os.path.join(root_folder, folder_path)
            
            if os.path.exists(absolute_path):
                count = 0
                for dirpath, dirnames, filenames in os.walk(absolute_path):
                    count += len(filenames)
                f.write(f"  {step_name}: {count} files\n")
            else:
                f.write(f"  {step_name}: FOLDER NOT FOUND\n")
    
    print("-" * 60)
    print("FILES SAVED:")
    print(f"  - {output_file} (simple list of unique filenames)")
    print(f"  - {detailed_output_file} (detailed CSV with full paths)")
    print(f"  - {summary_file} (scan summary)")
    print()
    print(f"Total files found: {total_files}")
    print(f"Unique filenames: {len(set(all_files))}")
    print(f"Folders processed: {folders_processed}/6")
    
    return all_files, detailed_files

# Alternative: Just list filenames without full recursion (top-level only)
def list_top_level_filenames():
    """
    List only top-level filenames (not going into subdirectories)
    """
    step_input_folder_name = []
    step_input_folder_name.append("Step1")
    step_input_folder_name.append("Step2")
    step_input_folder_name.append("Step3")
    step_input_folder_name.append("Step4_NEW")
    step_input_folder_name.append("Step5_New")
    step_input_folder_name.append("Step6_New")
    
    root_folder = os.getcwd()
    output_file = os.path.join(root_folder, "top_level_filenames.txt")
    
    print(f"Script running from: {root_folder}")
    print(f"Output will be saved to: {output_file}")
    print("-" * 60)
    
    all_files = []
    
    for step_name in step_input_folder_name:
        folder_path = os.path.join("data", "cluster", "endsim", step_name)
        absolute_path = os.path.join(root_folder, folder_path)
        
        if os.path.exists(absolute_path):
            files = [f for f in os.listdir(absolute_path) if os.path.isfile(os.path.join(absolute_path, f))]
            all_files.extend(files)
            print(f"{step_name}: {len(files)} files")
        else:
            print(f"{step_name}: Folder not found")
    
    with open(output_file, 'w') as f:
        for filename in sorted(set(all_files)):
            f.write(f"{filename}\n")
    
    print("-" * 60)
    print(f"Saved {len(set(all_files))} unique filenames to {output_file}")

if __name__ == "__main__":
    # Run the main recursive version
    print("Running recursive scan (includes subdirectories)...")
    print("=" * 60)
    all_files, detailed_files = list_filenames_in_folders()
    
    # If you only want top-level files, uncomment the line below:
    # print("\n" + "="*60)
    # print("Running top-level scan...")
    # list_top_level_filenames()