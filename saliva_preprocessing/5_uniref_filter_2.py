import aiohttp
import asyncio
import pandas as pd
import time


# use uniref database

async def fetch_unirefs(lists, timeout=10):

    # Create a list to store results
    # results = []

    # Create aiohttp client session
    async with aiohttp.ClientSession() as session:
        # Create tasks for concurrent requests
        tasks = []
        for id, single_list in enumerate(lists):
            task = asyncio.create_task(fetch_single_uniref(session, single_list, timeout))
            tasks.append(task)

        # Wait for all tasks to complete
        responses = await asyncio.gather(*tasks)

    # Process responses
    combined_result = sum(responses, [])

    return combined_result


async def fetch_single_uniref(session, mylist, timeout):

    queries = "%20OR%20".join(mylist)
    url = f"https://rest.uniprot.org/uniref/search?query={queries}"
    # Send GET request
    async with session.get(url, timeout=timeout) as response:
        # Check if request was successful
        if response.status == 200:
            # Parse JSON response
            data = await response.json()
            result = data["results"]
            output = []
            if len(result) > 0:
                for item in result:
                    output.append(item['representativeMember']["memberId"])
                    if "uniparcId" in item['representativeMember'].keys():
                        output.append(item['representativeMember']["uniparcId"])

                set1 = set(mylist)
                set2 = set(output)
                output = set1.intersection(set2)

            return list(output)
        else:
            print(f"Error fetching {url}: Status {response.status}")
            return None



def split_list_into_chunks(long_list, chunk_size=25):
    # Use list comprehension to break long_list into chunks of size chunk_size
    return [long_list[i:i + chunk_size] for i in range(0, len(long_list), chunk_size)]




# Example usage
async def main(input_df, output_file):

    uids = input_df['GeneShort'].tolist()
    uid_lists = split_list_into_chunks(uids, chunk_size=20)
    # Measure execution time
    start_time = time.time()
    # Run async API fetch
    results = await fetch_unirefs(uid_lists, timeout=600)
    # counts = [len(result) for result in results]
    # print(counts)
    results_df = pd.DataFrame({'ID': results})
    # Print results and execution time
    results_df.to_csv(output_file, index=False)
    print(f"\nExecution Time: {time.time() - start_time:.2f} seconds")




# Run the async main function
if __name__ == '__main__':

    unirefs90_plaque_df_2 = pd.read_csv("unirefs/uniref90_part2_plaque.csv")
    asyncio.run(main(unirefs90_plaque_df_2, "unirefs/output/uniref90_plaque_part2_subset.csv"))

    unirefs90_saliva_df_2 = pd.read_csv("unirefs/uniref90_part2_saliva.csv")
    asyncio.run(main(unirefs90_saliva_df_2, "unirefs/output/uniref90_saliva_part2_subset.csv"))
