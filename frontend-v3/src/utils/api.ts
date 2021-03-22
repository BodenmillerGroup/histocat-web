import ky from "ky";
import { apiUrl } from "../env";
import { AppToaster } from "./toaster";

export class ApiManager {
  private static _token: string;
  private static _api: typeof ky;

  static init(token: string) {
    ApiManager._token = token;
    ApiManager._api = ky.extend({
      prefixUrl: apiUrl,
      throwHttpErrors: true,
      hooks: {
        beforeRequest: [
          async (request) => {
            request.headers.set("Authorization", `Bearer ${ApiManager._token}`);
          },
        ],
        afterResponse: [
          async (request, options, response) => {
            if (!response.ok) {
              try {
                const errorJson = await response.json();
                return new Response(null, {
                  status: response.status,
                  statusText: errorJson.detail ? errorJson.detail : errorJson.error,
                });
              } catch (e) {
                console.log(e);
              }
            }
          },
        ],
      },
    });
  }

  static get api() {
    return ApiManager._api;
  }
}

export function displayApiError(error: any) {
  AppToaster.show({ message: error.message, intent: "danger" });
  if (error.response) {
    console.error("API error: ", error.response);
    if (error.response.status === 401) {
      // await this.actions.logOut();
    }
  }
}