<template>
  <v-col>
    <v-toolbar dense class="toolbar">
      <v-toolbar-title>Models</v-toolbar-title>
      <v-spacer />
      <v-toolbar-items>
        <v-btn v-if="isGroupAdmin" text :to="`/main/groups/${activeGroupId}/models/create`" color="primary">
          Add Model
        </v-btn>
      </v-toolbar-items>
    </v-toolbar>
    <v-card>
      <v-card-title>
        <v-spacer />
        <v-text-field v-model="search" append-icon="mdi-magnify" label="Search" single-line hide-details clearable />
      </v-card-title>
      <v-data-table
        :headers="headers"
        :items="items"
        :loading="!items"
        :search="search"
        :items-per-page="15"
        :footer-props="{
          itemsPerPageOptions: [10, 15, 20, -1],
          showFirstLastPage: true,
          showCurrentPage: true,
        }"
        multi-sort
      >
        <template v-slot:item.action="{ item }">
          <v-menu bottom left>
            <template v-slot:activator="{ on }">
              <v-btn icon v-on="on">
                <v-icon>mdi-dots-vertical</v-icon>
              </v-btn>
            </template>
            <v-list dense>
              <v-list-item
                :to="{
                  name: 'group-models-edit',
                  params: {
                    groupId: activeGroupId,
                    id: item.id,
                  },
                }"
              >
                <v-list-item-icon>
                  <v-icon color="grey">mdi-pencil-outline</v-icon>
                </v-list-item-icon>
                <v-list-item-content>
                  <v-list-item-title>Edit</v-list-item-title>
                </v-list-item-content>
              </v-list-item>
              <v-list-item v-if="isGroupAdmin" @click="deleteModel(item.id)">
                <v-list-item-icon>
                  <v-icon color="red accent-1">mdi-delete-outline</v-icon>
                </v-list-item-icon>
                <v-list-item-content>
                  <v-list-item-title>Delete</v-list-item-title>
                </v-list-item-content>
              </v-list-item>
            </v-list>
          </v-menu>
        </template>
      </v-data-table>
    </v-card>
  </v-col>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import { groupModule } from "@/modules/group";
import { modelsModule } from "@/modules/models";

@Component
export default class ModelsListView extends Vue {
  readonly groupContext = groupModule.context(this.$store);
  readonly modelsContext = modelsModule.context(this.$store);

  readonly headers = [
    {
      text: "Name",
      value: "name",
    },
    {
      text: "Description",
      value: "description",
    },
    {
      text: "Date",
      value: "created_at",
      filterable: false,
      width: "240",
    },
    {
      text: "Actions",
      value: "action",
      sortable: false,
      filterable: false,
      width: "70",
    },
  ];

  search = "";

  get activeGroupId() {
    return this.groupContext.getters.activeGroupId;
  }

  get isGroupAdmin() {
    return this.groupContext.getters.isGroupAdmin;
  }

  get items() {
    return this.modelsContext.getters.models;
  }

  async mounted() {
    await this.modelsContext.actions.getModels(+this.$router.currentRoute.params.groupId);
  }

  async deleteModel(id: number) {
    if (self.confirm("Are you sure you want to delete the model?")) {
      await this.modelsContext.actions.deleteModel(id);
    }
  }
}
</script>

<style scoped>
.toolbar {
  margin-bottom: 10px;
}
</style>
