import ky from "ky";

async function imageUrlToBase64(token: string, url: string): Promise<string> {
  const response = await ky.get(url, {
    headers: {
      Authorization: `Bearer ${token}`
    },
    timeout: false
  });
  const blob = await response.blob();
  const reader = new FileReader();
  return new Promise((resolve, reject) => {
    reader.onloadend = () => resolve(reader.result as string);
    reader.readAsDataURL(blob);
  });
}
