import os
import re
import pandas as pd
import numpy as np
from pathlib import Path

def parse_filename_direct(filename):
    """
    Parse all parameters from filename including silence parameters.
    Returns a dictionary with all parameters found.
    """
    base_name = os.path.splitext(filename)[0]
    
    # Initialize result dictionary
    params = {}
    
    # Extract basic parameters (always present)
    # epsM and epsSD
    eps_match = re.search(r'epsM([-\d.]+)_epsSD([-\d.]+)', base_name)
    if eps_match:
        params['epsM'] = float(eps_match.group(1))
        params['epsSD'] = float(eps_match.group(2))
    
    # Opinion parameters
    op_match = re.search(r'OpD([a-zA-Z-]+)_OpM([-\d.]+)_OpSD([-\d.]+)', base_name)
    if op_match:
        params['OpD'] = op_match.group(1)
        params['OpM'] = float(op_match.group(2))
        params['OpSD'] = float(op_match.group(3))
    
    # Network parameters
    net_match = re.search(r'Net([a-zA-Z-]+)___NAgents(\d+)', base_name)
    if net_match:
        params['NetworkType'] = net_match.group(1)
        params['NAgents'] = int(net_match.group(2))
    
    # Random seed
    rs_match = re.search(r'___RS(\d+)', base_name)
    if rs_match:
        params['RandomSeed'] = int(rs_match.group(1))
    
    # Media parameters
    media_match = re.search(r'__MedInfF([\d.]+)___MedD([a-zA-Z-]+)_MedN(\d+)_MedM([-\d.]+)_MedSD([-\d.]+)', base_name)
    if media_match:
        params['MedInfF'] = float(media_match.group(1))
        params['MedD'] = media_match.group(2)
        params['MedN'] = int(media_match.group(3))
        params['MedM'] = float(media_match.group(4))
        params['MedSD'] = float(media_match.group(5))
    
    # Silence parameters (optional - only in later steps)
    silence_match = re.search(r'_Silence_Alpha([-\d.]+)_Silence_Tau([-\d.]+)_Silence_Delta0([-\d.]+)_SilenceByBoundary\?([a-zA-Z]+)', base_name)
    if silence_match:
        params['Silence_Alpha'] = float(silence_match.group(1))
        params['Silence_Tau'] = float(silence_match.group(2))
        params['Silence_Delta0'] = float(silence_match.group(3))
        params['SilenceByBoundary'] = silence_match.group(4).lower() == 'true'
    else:
        # Set silence parameters as None for steps that don't have them
        params['Silence_Alpha'] = None
        params['Silence_Tau'] = None
        params['Silence_Delta0'] = None
        params['SilenceByBoundary'] = None
    
    return params

def get_all_filenames_from_folder(folder_path):
    """
    Get all CSV filenames from a folder
    """
    if not os.path.exists(folder_path):
        print(f"Warning: Folder does not exist: {folder_path}")
        return []
    
    filenames = [f for f in os.listdir(folder_path) if f.endswith('.csv')]
    return filenames

def determine_actual_parameters(df, min_unique_values=2):
    """
    Determine which parameters are ACTUALLY varied in the data.
    A parameter is considered a 'real variable' if it has at least 'min_unique_values' unique values.
    Parameters with only 1 unique value are constants and will be dropped.
    """
    param_columns = ['epsM', 'epsSD', 'OpM', 'OpSD', 'MedM', 'MedSD', 'MedInfF', 
                     'Silence_Alpha', 'Silence_Tau', 'Silence_Delta0']
    
    actual_params = {}
    
    for param in param_columns:
        if param in df.columns:
            unique_values = df[param].nunique()
            # Check if parameter exists (not all None) and has enough unique values
            has_values = df[param].notna().any()
            
            if has_values and unique_values >= min_unique_values:
                actual_params[param] = {
                    'unique_values': unique_values,
                    'values': df[param].unique(),
                    'min': df[param].min(),
                    'max': df[param].max()
                }
                print(f"  ✓ {param}: {unique_values} unique values (range: {df[param].min():.3f} - {df[param].max():.3f})")
            elif has_values and unique_values == 1:
                constant_value = df[param].iloc[0]
                print(f"  - {param}: CONSTANT ({constant_value}) - will be dropped")
            else:
                print(f"  ✗ {param}: NOT PRESENT or all NULL")
    
    return actual_params

def analyze_step_parameters(step_number, step_name, base_path="data/cluster/endsim"):
    """
    Analyze parameters for a single step by reading filenames directly
    """
    folder_path = os.path.join(base_path, step_name)
    print(f"\n{'='*70}")
    print(f"STEP {step_number}: {step_name}")
    print(f"Folder: {folder_path}")
    print(f"{'='*70}")
    
    # Get all filenames
    filenames = get_all_filenames_from_folder(folder_path)
    print(f"Found {len(filenames)} CSV files")
    
    if not filenames:
        return None
    
    # Parse all filenames
    parsed_data = []
    for filename in filenames:
        params = parse_filename_direct(filename)
        parsed_data.append(params)
    
    # Convert to DataFrame
    df = pd.DataFrame(parsed_data)
    
    # Display basic stats
    print(f"\n--- Parameter Analysis for Step {step_number} ---")
    
    # Determine which parameters are actually varied
    actual_params = determine_actual_parameters(df)
    
    # For silence parameters specifically, check if they exist
    print(f"\n--- Silence Parameter Status ---")
    if 'Silence_Tau' in df.columns:
        unique_tau = df['Silence_Tau'].nunique()
        non_null_tau = df['Silence_Tau'].notna().sum()
        print(f"  Silence_Tau: {unique_tau} unique values, {non_null_tau} non-null out of {len(df)}")
        if unique_tau > 1:
            print(f"    Values: {sorted(df['Silence_Tau'].unique())}")
        elif unique_tau == 1:
            print(f"    Constant value: {df['Silence_Tau'].iloc[0]}")
        else:
            print(f"    All values are NaN/None")
    else:
        print(f"  Silence_Tau column not found in parsed data")
    
    return df, actual_params

def analyze_all_steps():
    """
    Analyze parameters for all steps (1-6)
    """
    step_input_folder_name = [
        "Step1", "Step2", "Step3", 
        "Step4_NEW", "Step5_New", "Step6_New"
    ]
    
    step_titles = [
        "Simple HK with Full Network",
        "Simple HK with Scale-Free Network",
        "Heterogenous Openness",
        "Heterogenous Openness with normal opinion distributions",
        "Media and Influence",
        "Media with Silence"
    ]
    
    results = {}
    
    for step_no, (step_name, step_title) in enumerate(zip(step_input_folder_name, step_titles), start=1):
        print(f"\n{'#'*70}")
        print(f"# Analyzing Step {step_no}: {step_title}")
        print(f"# Folder: {step_name}")
        print(f"{'#'*70}")
        
        result = analyze_step_parameters(step_no, step_name)
        
        if result:
            df, actual_params = result
            results[step_no] = {
                'step_name': step_name,
                'step_title': step_title,
                'df': df,
                'actual_params': actual_params,
                'n_simulations': len(df)
            }
    
    # Summary across all steps
    print(f"\n\n{'='*70}")
    print("SUMMARY ACROSS ALL STEPS")
    print(f"{'='*70}")
    
    summary_data = []
    for step_no in results:
        params_with_variation = list(results[step_no]['actual_params'].keys())
        has_silence = any('Silence' in p for p in params_with_variation)
        
        summary_data.append({
            'Step': step_no,
            'Step Name': results[step_no]['step_title'],
            'N_Simulations': results[step_no]['n_simulations'],
            'Parameters_With_Variation': len(params_with_variation),
            'Has_Silence_Parameters': has_silence,
            'Silence_Tau_Values': results[step_no]['df']['Silence_Tau'].unique() if 'Silence_Tau' in results[step_no]['df'].columns else []
        })
        
        print(f"\nStep {step_no}: {results[step_no]['step_title']}")
        print(f"  Simulations: {results[step_no]['n_simulations']}")
        print(f"  Parameters with variation: {params_with_variation}")
        if has_silence:
            tau_vals = results[step_no]['df']['Silence_Tau'].unique()
            print(f"  Silence_Tau values: {sorted([v for v in tau_vals if pd.notna(v)])}")
    
    # Create summary DataFrame
    summary_df = pd.DataFrame(summary_data)
    print(f"\n{'='*70}")
    print("SUMMARY TABLE:")
    print(summary_df.to_string(index=False))
    
    return results, summary_df

def verify_silence_parsing():
    """
    Test function to verify silence parameter parsing with example filename
    """
    test_filename = "EndSim_epsM0.21_epsSD0.05___OpDnormal_OpM0.25_OpSD0___NetScale-free___NAgents1000___RS2__MedInfF1___MedDdeterministic-normal_MedN10_MedM-0.75_MedSD0.7_Silence_Alpha0.8_Silence_Tau1_Silence_Delta00.75_SilenceByBoundary?true.csv"
    
    print("Testing silence parameter parsing...")
    params = parse_filename_direct(test_filename)
    
    print(f"\nParsed parameters:")
    for key, value in params.items():
        if 'Silence' in key:
            print(f"  {key}: {value} (type: {type(value)})")
        else:
            print(f"  {key}: {value}")
    
    # Verify silence values
    expected_silence = {
        'Silence_Alpha': 0.8,
        'Silence_Tau': 1.0,
        'Silence_Delta0': 0.75,
        'SilenceByBoundary': True
    }
    
    print(f"\nVerification:")
    all_correct = True
    for key, expected in expected_silence.items():
        actual = params.get(key)
        if actual == expected:
            print(f"  ✓ {key}: {actual} (correct)")
        else:
            print(f"  ✗ {key}: expected {expected}, got {actual}")
            all_correct = False
    
    if all_correct:
        print("\n✅ Silence parameters parsed correctly!")
    else:
        print("\n❌ Silence parameters not parsed correctly")
    
    return params

if __name__ == "__main__":
    # First test the silence parsing with example filename
    print("Testing silence parameter parsing with example filename...")
    verify_silence_parsing()
    
    print("\n\n")
    print("="*70)
    print("ANALYZING ALL STEPS")
    print("="*70)
    
    # Analyze all steps
    results, summary_df = analyze_all_steps()
    
    # Save results to file
    output_file = "parameter_analysis_summary.csv"
    summary_df.to_csv(output_file, index=False)
    print(f"\n✅ Summary saved to {output_file}")
    
    # Specifically check Step 6 silence parameters
    if 6 in results:
        print(f"\n{'='*70}")
        print("STEP 6 SILENCE PARAMETER DETAILS")
        print(f"{'='*70}")
        step6_df = results[6]['df']
        
        print(f"\nSilence_Tau distribution in Step 6:")
        tau_counts = step6_df['Silence_Tau'].value_counts()
        for tau_value, count in tau_counts.items():
            print(f"  Silence_Tau = {tau_value}: {count} simulations ({count/len(step6_df)*100:.1f}%)")
        
        print(f"\nSilence_Alpha distribution:")
        alpha_counts = step6_df['Silence_Alpha'].value_counts()
        for alpha_value, count in alpha_counts.items():
            print(f"  Silence_Alpha = {alpha_value}: {count} simulations")
        
        print(f"\nSilence_Delta0 distribution:")
        delta0_counts = step6_df['Silence_Delta0'].value_counts()
        for delta0_value, count in delta0_counts.items():
            print(f"  Silence_Delta0 = {delta0_value}: {count} simulations")