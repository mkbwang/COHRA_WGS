import aiohttp
import asyncio
import pandas as pd
import time
import numpy as np


# use uniprotkb database

async def fetch_unirefs(lists, timeout=10):

    # Create a list to store results
    # results = []

    # Create aiohttp client session
    async with aiohttp.ClientSession() as session:
        # Create tasks for concurrent requests
        tasks = []
        for id, single_list in enumerate(lists):
            queries = "%20OR%20".join(single_list)
            url = f"https://rest.uniprot.org/uniprotkb/search?query={queries}"
            task = asyncio.create_task(fetch_single_uniref(session, url, timeout, id))
            tasks.append(task)

        # Wait for all tasks to complete
        responses = await asyncio.gather(*tasks)

        # Process responses
    combined_result = pd.concat(responses, ignore_index=True, axis=0)

    return combined_result


async def fetch_single_uniref(session, url, timeout, id):

    # Send GET request
    async with session.get(url, timeout=timeout) as response:
        # Check if request was successful
        if response.status == 200:
            # Parse JSON response
            data = await response.json()
            result = data["results"]
            existence = [item["entryType"] for item in result]
            ID = [item["primaryAccession"] for item in result]
            existence_df = pd.DataFrame({"ID": ID, "Existence": existence})
            return existence_df
        else:
            print(f"Error fetching {url}: Status {response.status}")
            return None



def split_list_into_chunks(long_list, chunk_size=25):
    # Use list comprehension to break long_list into chunks of size chunk_size
    return [long_list[i:i + chunk_size] for i in range(0, len(long_list), chunk_size)]




# Example usage
async def main(input_df, output_file):

    uids = input_df['GeneShort'].tolist()
    uid_lists = split_list_into_chunks(uids)
    # Measure execution time
    start_time = time.time()
    # Run async API fetch
    results = await fetch_unirefs(uid_lists, timeout=600)
    # Print results and execution time
    results.to_csv(output_file, index=False)
    print(f"\nExecution Time: {time.time() - start_time:.2f} seconds")




# Run the async main function
if __name__ == '__main__':

    print("processing plaque...")
    unirefs90_plaque_df_1 = pd.read_csv("unirefs/uniref90_part1_plaque.csv")
    num_samples = len(unirefs90_plaque_df_1)
    samples_lists = split_list_into_chunks(np.arange(num_samples), chunk_size=150000)
    for i, sample_list in enumerate(samples_lists):
        print(i)
        unirefs90_plaque_df_subset = unirefs90_plaque_df_1.iloc[sample_list, :]
        asyncio.run(main(unirefs90_plaque_df_subset, f"unirefs/output/uniref90_plaque_part1_subset_{i}.csv"))

    print("processing saliva...")
    unirefs90_saliva_df_1 = pd.read_csv("unirefs/uniref90_part1_saliva.csv")
    num_samples = len(unirefs90_saliva_df_1)
    samples_lists = split_list_into_chunks(np.arange(num_samples), chunk_size=150000)
    for i, sample_list in enumerate(samples_lists):
        print(i)
        unirefs90_saliva_df_subset = unirefs90_saliva_df_1.iloc[sample_list, :]
        asyncio.run(main(unirefs90_saliva_df_subset, f"unirefs/output/uniref90_saliva_part1_subset_{i}.csv"))
