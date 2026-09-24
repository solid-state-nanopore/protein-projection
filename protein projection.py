import argparse
import math
import numpy as np

def read_ca_coordinates(pdb_file):
    coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                atom_name = line[12:16].strip()
                if atom_name == 'CA':
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append([x, y, z])
                    except:
                        continue
    return np.array(coords)

def calculate_center(coords):
    return np.mean(coords, axis=0)

def rotate_3d(coords, angle_x, angle_y, center):
    homogenous_coords = np.hstack([coords - center, np.ones((len(coords), 1))])
    
    theta_x = math.radians(angle_x)
    rot_x = np.array([
        [1, 0, 0, 0],
        [0, math.cos(theta_x), -math.sin(theta_x), 0],
        [0, math.sin(theta_x), math.cos(theta_x), 0],
        [0, 0, 0, 1]
    ])
    
    theta_y = math.radians(angle_y)
    rot_y = np.array([
        [math.cos(theta_y), 0, math.sin(theta_y), 0],
        [0, 1, 0, 0],
        [-math.sin(theta_y), 0, math.cos(theta_y), 0],
        [0, 0, 0, 1]
    ])
    
    transformed = homogenous_coords @ rot_x.T @ rot_y.T
    
    return transformed[:, :3] + center

def calculate_projected_area(proj_2d, step_size, atom_radius=1.7):
    if len(proj_2d) == 0:
        return 0.0
    
    min_x = np.min(proj_2d[:, 0]) - atom_radius
    max_x = np.max(proj_2d[:, 0]) + atom_radius
    min_y = np.min(proj_2d[:, 1]) - atom_radius
    max_y = np.max(proj_2d[:, 1]) + atom_radius
    
    x_grid = np.arange(min_x, max_x + step_size, step_size)
    y_grid = np.arange(min_y, max_y + step_size, step_size)
    grid = np.zeros((len(x_grid), len(y_grid)), dtype=bool)
    
    for x, y in proj_2d:
        x_indices = np.where((x_grid >= x - atom_radius) & (x_grid <= x + atom_radius))[0]
        y_indices = np.where((y_grid >= y - atom_radius) & (y_grid <= y + atom_radius))[0]
        
        xx, yy = np.meshgrid(x_grid[x_indices], y_grid[y_indices], indexing='ij')
        
        dist_sq = (xx - x)**2 + (yy - y)**2
        covered = dist_sq <= atom_radius**2
        
        grid[np.ix_(x_indices, y_indices)] |= covered
    
    return np.sum(grid) * (step_size ** 2)

def main():
    parser = argparse.ArgumentParser(description='Protein dual-axis rotation projected area analysis')
    parser.add_argument('input_pdb', help='Input PDB file path')
    parser.add_argument('-o', '--output', default='area_matrix.dat',
                      help='Output filename (default: area_matrix.dat)')
    parser.add_argument('--step', type=float, default=1.0,
                      help='Grid step size (Å, default: 1.0)')
    parser.add_argument('--radius', type=float, default=1.7,
                      help='Atom radius (Å, default: 1.7)')
    
    args = parser.parse_args()
    
    coords = read_ca_coordinates(args.input_pdb)
    if len(coords) == 0:
        print("Error: No CA atoms found")
        return
    center = calculate_center(coords)
    print(f"Center coordinates: {center}")
    
    matrix = np.zeros((360, 360), dtype=np.float32)
    
    for angle_x in np.arange(0, 180, 0.5):
        for angle_y in np.arange(0, 180, 0.5):
            angle_x_index = int(angle_x * 2)
            angle_y_index = int(angle_y * 2)
            rotated = rotate_3d(coords, angle_x, angle_y, center)
            
            proj_2d = rotated[:, :2]
            
            area = calculate_projected_area(proj_2d, args.step, args.radius)
            matrix[angle_x_index, angle_y_index] = area
            
            print(f"X:{angle_x:3.1f}° Y:{angle_y:3.1f}° => Area: {area:.2f} Å²")

    
    np.savetxt(args.output, matrix, 
              fmt='%.4f',
              delimiter='\t',
              header=' '.join([f'Y{deg:03.1f}' for deg in np.arange(0, 180, 0.5)]))
    
    print(f"Results saved to {args.output}")

if __name__ == '__main__':
    main()
