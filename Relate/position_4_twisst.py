import sys

# Function to generate intervals between numbers
def generate_intervals(file_path):
    with open(file_path, 'r') as file:
        numbers = [int(line.strip()) for line in file]

    for i in range(len(numbers) - 1):
        start = numbers[i]
        end = numbers[i + 1] - 1
        print(f"{start} {end}")

    # For the last number, add 1000
    last_number = numbers[-1]
    print(f"{last_number} {last_number + 1000}")

# Main part
if __name__ == "__main__":
    file_path = sys.argv[1]
    generate_intervals(file_path)
