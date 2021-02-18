import { WebSocketMessage } from "./WebSocketMessage";
import { AppToaster } from "./toaster";
import { useProjectsStore } from "../modules/projects";
import { useDatasetsStore } from "../modules/datasets";

export class WebSocketManager {
  static socket: WebSocket;
  static token: string;
  static protocol: string;

  static init(token: string) {
    WebSocketManager.token = token;
    WebSocketManager.protocol = window.location.protocol === "https:" ? "wss:" : "ws:";
  }

  static connect(projectId: number) {
    // Dispose previous connections due to login/logout
    WebSocketManager.close();
    const url = `${WebSocketManager.protocol}//${window.location.hostname}/ws/${projectId}?token=${WebSocketManager.token}`;

    WebSocketManager.socket = new WebSocket(url);
    WebSocketManager.socket.onopen = (event: Event) => {
      console.log(`WebSocket connection to project ${projectId} opened`);
    };

    WebSocketManager.socket.onclose = (event: CloseEvent) => {
      console.log("WebSocket is closed. Reconnecting...", event.reason);
      setTimeout(() => {
        WebSocketManager.connect(projectId);
      }, 1000);
    };

    WebSocketManager.socket.onmessage = async (event: MessageEvent) => {
      const json = JSON.parse(event.data);
      const message = new WebSocketMessage(json);
      console.log(message);
      const activeProjectId = useProjectsStore.getState().activeProjectId;
      if (message.projectId === activeProjectId) {
        switch (message.type) {
          case "slide_imported": {
            const getProjectData = useProjectsStore.getState().getProjectData;
            await getProjectData(message.projectId);
            AppToaster.show({ message: "Slide successfully imported", intent: "success" });
            break;
          }
          case "segmentation_ready":
          case "dataset_imported": {
            const getProjectDatasets = useDatasetsStore.getState().getProjectDatasets;
            await getProjectDatasets(message.projectId);
            AppToaster.show({ message: "Dataset successfully imported", intent: "success" });
            break;
          }
          case "result_ready": {
            const activeDatasetId = useDatasetsStore.getState().activeDatasetId;
            if (message.payload.dataset_id === activeDatasetId) {
              // WebSocketManager.resultContext.mutations.addEntity(message.payload);
              AppToaster.show({ message: "Pipeline processing result is ready", intent: "success" });
            }
            break;
          }
        }
      }
    };

    WebSocketManager.socket.onerror = (event: Event) => {
      console.log("WebSocket error: ", event);
    };
  }

  static close() {
    if (WebSocketManager.socket) {
      WebSocketManager.socket.onopen = null;
      WebSocketManager.socket.onclose = null;
      WebSocketManager.socket.onmessage = null;
      WebSocketManager.socket.onerror = null;
      WebSocketManager.socket.close();
      console.log("WebSocket connection closed");
    }
  }
}
