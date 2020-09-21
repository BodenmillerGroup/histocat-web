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
import { resultsModule, ResultsState } from "@/modules/results";
import { ResultsGetters } from "@/modules/results/getters";
import { ResultsMutations } from "@/modules/results/mutations";
import { ResultsActions } from "@/modules/results/actions";

export class WebSocketManager {
  static mainContext: Context<Module<MainState, MainGetters, MainMutations, MainActions>>;
  static projectsContext: Context<Module<ProjectsState, ProjectsGetters, ProjectsMutations, ProjectsActions>>;
  static datasetContext: Context<Module<DatasetsState, DatasetsGetters, DatasetsMutations, DatasetsActions>>;
  static resultContext: Context<Module<ResultsState, ResultsGetters, ResultsMutations, ResultsActions>>;
  static socket: WebSocket;
  static token: string;
  static protocol: string;

  static init(store: Store<any>) {
    WebSocketManager.mainContext = mainModule.context(store);
    WebSocketManager.projectsContext = projectsModule.context(store);
    WebSocketManager.datasetContext = datasetsModule.context(store);
    WebSocketManager.resultContext = resultsModule.context(store);
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
          case "dataset_imported": {
            WebSocketManager.datasetContext.actions.getProjectDatasets(message.projectId);
            WebSocketManager.mainContext.mutations.addNotification({
              content: "Dataset successfully imported",
              color: "success",
            });
            break;
          }
          case "tsne_result_ready": {
            WebSocketManager.resultContext.mutations.addEntity(message.payload);
            WebSocketManager.mainContext.mutations.addNotification({
              content: "t-SNE result is ready",
              color: "success",
            });
            break;
          }
          case "umap_result_ready": {
            WebSocketManager.resultContext.mutations.addEntity(message.payload);
            WebSocketManager.mainContext.mutations.addNotification({
              content: "UMAP result is ready",
              color: "success",
            });
            break;
          }
          case "phenograph_result_ready": {
            WebSocketManager.resultContext.mutations.addEntity(message.payload);
            WebSocketManager.mainContext.mutations.addNotification({
              content: "PhenoGraph result is ready",
              color: "success",
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
