const envApiUrl = `${process.env.VUE_APP_PROTOCOL}://${process.env.VUE_APP_DOMAIN}`;

export const apiUrl = `${envApiUrl}/api/v1`;
export const appName = process.env.VUE_APP_NAME;
