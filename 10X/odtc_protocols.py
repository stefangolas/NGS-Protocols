import xml.etree.ElementTree as ET
from datetime import datetime
from pyhamilton.odtc import ThermalCyclerProtocol

def calculate_pcr_cycles(cdna_mass_ng):
    """Calculate optimal PCR cycles based on cDNA input mass."""
    # (max_mass, cycles) tuples in ascending order
    mass_cycle_ranges = [
        (50, 15),     # 0.25-50 ng → 14-16 cycles (use 15)
        (250, 13),    # 50-250 ng → 12-14 cycles (use 13)  
        (600, 11),    # 250-600 ng → 10-12 cycles (use 11)
        (1100, 9),    # 600-1,100 ng → 8-10 cycles (use 9)
        (1500, 7),    # 1,100-1,500 ng → 6-8 cycles (use 7)
        (float('inf'), 5)  # >1,500 ng → 5 cycles
    ]
    
    for max_mass, cycles in mass_cycle_ranges:
        if cdna_mass_ng <= max_mass:
            return cycles


def create_fragmentation_protocol():
    """
    Protocol 3.1: GEX Fragmentation, End Repair & A-tailing
    - Pre-cool block to 4°C
    - Fragmentation: 32°C for 5 min
    - End Repair & A-Tailing: 65°C for 30 min
    - Hold at 4°C
    """
    protocol = ThermalCyclerProtocol(creator="10x_genomics", fluid_quantity=50)
    protocol.set_pre_method(block_temp=4, lid_temp=65)
    
    # Fragmentation step: 32°C for 5 minutes (300 seconds)
    protocol.add_step(plateau_temp=32, plateau_time=300, lid_temp=65)
    
    # End Repair & A-Tailing: 65°C for 30 minutes (1800 seconds)
    protocol.add_step(plateau_temp=65, plateau_time=1800, lid_temp=65)
    
    # Final hold at 4°C
    protocol.add_step(plateau_temp=4, plateau_time=600, lid_temp=65)
    
    return protocol


def create_adapter_ligation_protocol():
    """
    Protocol 3.3: GEX Adaptor Ligation
    - 20°C for 15 min
    - Hold at 4°C
    """
    protocol = ThermalCyclerProtocol(creator="10x_genomics", fluid_quantity=100)
    protocol.set_pre_method(block_temp=20, lid_temp=30)
    
    # Ligation step: 20°C for 15 minutes (900 seconds)
    protocol.add_step(plateau_temp=20, plateau_time=900, lid_temp=30)
    
    # Hold at 4°C
    protocol.add_step(plateau_temp=4, plateau_time=600, lid_temp=30)
    
    return protocol


def create_sample_index_pcr_protocol(cdna_mass_ng=None, num_cycles=None):
    """
    Protocol 3.5: GEX Sample Index PCR
    - Initial denaturation: 98°C for 45 sec
    - Cycling (variable cycles based on cDNA input):
      - 98°C for 20 sec
      - 54°C for 30 sec  
      - 72°C for 20 sec
    - Final extension: 72°C for 1 min
    - Hold at 4°C
    
    Args:
        cdna_mass_ng: cDNA input mass in ng (auto-calculates cycles)
        num_cycles: Manual cycle number override
    """
    protocol = ThermalCyclerProtocol(creator="10x_genomics", fluid_quantity=100)
    protocol.set_pre_method(lid_temp=105)
    
    # Determine cycle number
    if num_cycles is not None:
        cycles = num_cycles
    elif cdna_mass_ng is not None:
        cycles = calculate_pcr_cycles(cdna_mass_ng)
    else:
        cycles = 12  # Default for ~100 ng input
    
    # Initial denaturation: 98°C for 45 seconds
    protocol.add_step(plateau_temp=98, plateau_time=45, lid_temp=105)
    
    # PCR cycling using existing add_pcr_cycle method
    protocol.add_pcr_cycle(
        denaturation_temp=98, denaturation_time=20,
        annealing_temp=54, annealing_time=30,
        extension_temp=72, extension_time=20,
        num_cycles=cycles
    )
    
    # Final extension using existing add_final_extension method
    protocol.add_final_extension(temp=72, time=60)
    
    # Hold at 4°C
    protocol.add_step(plateau_temp=4, plateau_time=600, lid_temp=105)
    
    return protocol


def get_cdna_mass_from_user():
    """Prompt user for cDNA mass and validate input."""
    while True:
        try:
            mass_str = input("\nEnter your total cDNA yield (ng): ")
            mass = float(mass_str)
            if mass <= 0:
                print("Please enter a positive number.")
                continue
            return mass
        except ValueError:
            print("Please enter a valid number.")
        except KeyboardInterrupt:
            print("\nOperation cancelled.")
            return None




if __name__ == "__main__":
    print("=== 10x Genomics Thermal Cycler Protocol Generator ===")
    print("\nThis script generates thermal cycler protocols for:")
    print("1. Fragmentation, End Repair & A-tailing")
    print("2. Adapter Ligation")
    print("3. Sample Index PCR (optimized for your cDNA input)")
    
    # Get cDNA mass from user
    cdna_mass = get_cdna_mass_from_user()
    
    if cdna_mass is None:
        exit()
    
    # Calculate cycles
    pcr_cycles = calculate_pcr_cycles(cdna_mass)    
    print(f"\nGenerating protocols with {pcr_cycles} PCR cycles...")
    
    # Generate all protocols
    try:
        # Protocol 1: Fragmentation
        frag_protocol = create_fragmentation_protocol()
        frag_filename = "10x_fragmentation_protocol.xml"
        frag_protocol.generate_xml(frag_filename)
        print(f"✓ Generated: {frag_filename}")
        
        # Protocol 2: Adapter Ligation
        ligation_protocol = create_adapter_ligation_protocol()
        ligation_filename = "10x_adapter_ligation_protocol.xml"
        ligation_protocol.generate_xml(ligation_filename)
        print(f"✓ Generated: {ligation_filename}")
        
        # Protocol 3: Sample Index PCR
        pcr_protocol = create_sample_index_pcr_protocol(cdna_mass_ng=cdna_mass, num_cycles=pcr_cycles)
        pcr_filename = f"10x_sample_index_pcr.xml"
        pcr_protocol.generate_xml(pcr_filename)
        print(f"✓ Generated: {pcr_filename}")
        
        print(f"\n=== Protocol Generation Complete! ===")
        print(f"\nSummary:")
        print(f"  cDNA Input: {cdna_mass} ng")
        print(f"  PCR Cycles: {pcr_cycles}")
        print(f"  Files generated: 3")
        
        print(f"\nProtocol execution order:")
        print(f"  1. {frag_filename}")
        print(f"  2. {ligation_filename}")
        print(f"  3. {pcr_filename}")
        
    except Exception as e:
        print(f"\nError generating protocols: {e}")
        print("Please check that the ThermalCyclerProtocol class is properly imported.")