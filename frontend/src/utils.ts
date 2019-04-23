export const getLocalToken = () => localStorage.getItem('token');

export const saveLocalToken = (token: string) => localStorage.setItem('token', token);

export const removeLocalToken = () => localStorage.removeItem('token');

export function authHeaders(token: string, contentType?: string) {
  if (contentType) {
    return {
      headers: {
        Authorization: `Bearer ${token}`,
        'Content-Type': contentType
      },
    };
  } else {
    return {
      headers: {
        Authorization: `Bearer ${token}`,
      },
    };
  }
}
