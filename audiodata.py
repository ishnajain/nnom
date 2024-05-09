import wave
import csv
import struct

def wav_to_csv(input_wav_file, output_csv_file):
    # Open the WAV file
    with wave.open(input_wav_file, 'rb') as wav_file:
        # Get the parameters of the WAV file
        channels = wav_file.getnchannels()
        sample_width = wav_file.getsampwidth()
        frame_rate = wav_file.getframerate()
        num_frames = wav_file.getnframes()

        # Make sure the WAV file is mono with 16 kHz sample rate and 16-bit PCM
        if channels != 1 or sample_width != 2 or frame_rate != 16000:
            print("Unsupported WAV file format. Please ensure the WAV file is mono with 16 kHz sample rate and 16-bit PCM.")
            return

        # Read all frames from the WAV file
        frames = wav_file.readframes(num_frames)

        # Extract the raw PCM data as a list of 16-bit integers
        raw_pcm_data = list(struct.unpack("<" + "h" * (len(frames) // 2), frames))

        # Write the raw PCM data to a CSV file
        with open(output_csv_file, 'w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            # Write all samples as a single row in the CSV file
            csv_writer.writerow(raw_pcm_data)

# Usage example
input_wav_file = "right1.wav"
output_csv_file = "right1.csv"
wav_to_csv(input_wav_file, output_csv_file)
