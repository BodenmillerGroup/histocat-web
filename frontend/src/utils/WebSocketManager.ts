import { datasetsModule, DatasetsState } from "@/modules/datasets";
import { DatasetsActions } from "@/modules/datasets/actions";
import { DatasetsGetters } from "@/modules/datasets/getters";
import { DatasetsMutations } from "@/modules/datasets/mutations";
import { projectsModule, ProjectsState } from "@/modules/projects";
import { ProjectsActions } from "@/modules/projects/actions";
import { ProjectsGetters } from "@/modules/projects/getters";
import { ProjectsMutations } from "@/modules/projects/mutations";
import { mainModule, MainState } from "@/modules/main";
import { MainActions } from "@/modules/main/actions";
import { MainGetters } from "@/modules/main/getters";
import { MainMutations } from "@/modules/main/mutations";
import { WebSocketMessage } from "@/utils/WebSocketMessage";
import { Store } from "vuex";
import { Context, Module } from "vuex-smart-module";
import { cellsModule, CellsState } from "@/modules/cells";
import { CellsGetters } from "@/modules/cells/getters";
import { CellsMutations } from "@/modules/cells/mutations";
import { CellsActions } from "@/modules/cells/actions";

export class WebSocketManager {
  static mainContext: Context<Module<MainState, MainGetters, MainMutations, MainActions>>;
  static projectsContext: Context<Module<ProjectsState, ProjectsGetters, ProjectsMutations, ProjectsActions>>;
  static datasetContext: Context<Module<DatasetsState, DatasetsGetters, DatasetsMutations, DatasetsActions>>;
  static cellsContext: Context<Module<CellsState, CellsGetters, CellsMutations, CellsActions>>;
  static socket: WebSocket;
  static token: string;
  static protocol: string;

  static init(store: Store<any>) {
    WebSocketManager.mainContext = mainModule.context(store);
    WebSocketManager.projectsContext = projectsModule.context(store);
    WebSocketManager.datasetContext = datasetsModule.context(store);
    WebSocketManager.cellsContext = cellsModule.context(store);
    WebSocketManager.token = WebSocketManager.mainContext.getters.token;
    WebSocketManager.protocol = self.location.protocol === "https:" ? "wss:" : "ws:";
  }

  static connect(projectId: number) {
    // Dispose previous connections due to login/logout
    WebSocketManager.close();
    const url = `${WebSocketManager.protocol}//${location.hostname}/ws/${projectId}?token=${WebSocketManager.token}`;

    WebSocketManager.socket = new WebSocket(url);
    WebSocketManager.socket.onopen = (event: Event) => {
      console.log("WebSocket connection opened");
    };

    WebSocketManager.socket.onclose = (event: CloseEvent) => {
      console.log("WebSocket is closed. Reconnecting...", event.reason);
      setTimeout(() => {
        WebSocketManager.connect(projectId);
      }, 1000);
    };

    WebSocketManager.socket.onmessage = (event: MessageEvent) => {
      const json = JSON.parse(event.data);
      const message = new WebSocketMessage(json);
      console.log(message);
      if (message.projectId === WebSocketManager.projectsContext.getters.activeProjectId) {
        switch (message.type) {
          case "slide_imported": {
            WebSocketManager.projectsContext.actions.getProjectData(message.projectId);
            WebSocketManager.mainContext.mutations.addNotification({
              content: "Slide successfully imported",
              color: "success",
            });
            break;
          }
          case "segmentation_ready":
          case "dataset_imported": {
            WebSocketManager.datasetContext.actions.getProjectDatasets(message.projectId);
            WebSocketManager.mainContext.mutations.addNotification({
              content: "Dataset is ready",
              color: "success",
            });
            break;
          }
          case "result_ready": {
            if (message.payload.dataset_id === WebSocketManager.datasetContext.getters.activeDatasetId) {
              WebSocketManager.cellsContext.mutations.addEntity(message.payload);
              WebSocketManager.mainContext.mutations.addNotification({
                content: "Pipeline processing complete",
                color: "success",
              });
            }
            break;
          }
          case "error": {
            WebSocketManager.mainContext.mutations.addNotification({
              content: message.payload,
              color: "error",
            });
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
    }
  }
}
