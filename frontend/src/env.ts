const env = process.env.VUE_APP_ENV;

let envApiUrl = '';
let wsApiUrl = '';

if (env === 'production') {
  envApiUrl = `http://${process.env.VUE_APP_DOMAIN_PROD}`;
  wsApiUrl = `ws://${process.env.VUE_APP_DOMAIN_PROD}`;
} else if (env === 'staging') {
  envApiUrl = `http://${process.env.VUE_APP_DOMAIN_STAG}`;
  wsApiUrl = `ws://${process.env.VUE_APP_DOMAIN_STAG}`;
} else {
  envApiUrl = `http://${process.env.VUE_APP_DOMAIN_DEV}`;
  wsApiUrl = `ws://${process.env.VUE_APP_DOMAIN_DEV}`;
}

export const apiUrl = envApiUrl;
export const wsUrl = wsApiUrl;
export const appName = process.env.VUE_APP_NAME;
